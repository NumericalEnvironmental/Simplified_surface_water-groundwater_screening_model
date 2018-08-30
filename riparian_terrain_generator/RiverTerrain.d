
// **********************************************************************************************
//
// RiverTerrain.d - a terrain generator for an eroded river valley
//
// by Walt McNab
//
// Steps:
// (1) Posit river network shape in 3-D
// (2) Shape terrain surrounding drainages, based on stream proximity
//
// **********************************************************************************************


// the usual libraries ...
import std.stdio;           // I/O and file system
import std.math;            // math function library
import std.algorithm;       // algorithm for working with ranges, sequences, etc.
import std.array;           // array operations
import std.string;          // string function support
import std.conv;            // automatic conversions between types
import std.typecons;        // misc. meta-programming support, including type constructors, tuples, etc.
import std.range;           // enables range functions
import std.mathspecial;     // for normal and inverse normal distribution functions
import std.random;          // for random number generation


// misc. support functions


double CumNorm(double x, double avg, double stdev){
    // return the cumulative distribution percentile
    double xStandard = (x - avg)/stdev;
    return normalDistribution(xStandard);
    }


double BoundedRandUniform(double avg, double stdev, double valMin, double valMax){
    // random sample from normal distribution, bounded on both ends
    double p;
    double r = 1e+10;
    while((r < valMin) || (r > valMax)){
        p = uniform01();
        r = avg + stdev*normalDistributionInverse(p);
        }
    return r;
    }


double Distance(double x1, double y1, double x2, double y2){
    // distance between (x1, y1) and (x2, y2)
    return sqrt((x1-x2)^^2. + (y1-y2)^^2.);
    }    
    
    
void WriteReaches(Reach[] reach){
    // write reaches properties to output file
    string header, rchLine;
    auto outputFile = File("reaches.csv", "w"); 
    header = "x" ~ "," ~ "y" ~ "," ~ "z"  ~ "," ~ "L";
    header ~= ","~ "w"  ~ "," ~ "orient" ~ "," ~ "h";
    header ~= "," ~ "n" ~ "," ~ "K" ~ "," ~ "fixed";
    outputFile.writeln(header);
    foreach (rch; reach){
        rchLine = to!string(rch.x) ~ "," ~ to!string(rch.y) ~ "," ~ to!string(rch.z)  ~ "," ~ to!string(rch.L);
        rchLine ~= "," ~ to!string(rch.w)  ~ "," ~ to!string(rch.orient) ~ "," ~ to!string(rch.h);
        rchLine ~= "," ~ to!string(rch.n) ~ "," ~ to!string(rch.K) ~ "," ~ to!string(rch.fixed);
        outputFile.writeln(rchLine);        
        }
    outputFile.close();
    writeln("Wrote results to reaches file.");
    }
    
    
// classes  
    

class Point{
    double x, y, z;
    this(double x, double y, double z){
        this.x = x;
        this.y = y;
        this.z = z;
        }    
    }


class Reach: Point{

    // reach parameters (to be written separately to file for RiverRun.py)
    int rank;
    double L, w, orient, h, n, K;
    double zStart, zEnd;
    bool fixed;
    
    // constructor
    this(double x, double y, double z, double L, double w, double orient,
        double n, double K, bool fixed, int rank, double zStart, double zEnd){
        super(x, y, z);             // Reach class inherits location and elevation from Point class        
        this.L = L;                 // reach length (default; used to write file for RiverRun, etc.)
        this.w = w;                 // reach width 
        this.orient = orient;       // directional orientation, radians N of due E (i.e., counter-clockwise)
        this.n = n;                 // Manning coefficient (default; used to write file for RiverRun, etc.)
        this.K = K;                 // stream bed conductivity (default; used to write file for RiverRun, etc.)
        this.fixed = fixed;         // boolean flag indicating fixed-water level condition
        this.h = 0.;                // water depth (default; used to write file for RiverRun, etc.)
        this.rank = rank;           // no. or breaks off of main stream, as tributary
        this.zStart = zStart;       // starting and ending elevations
        this.zEnd = zEnd;
        }
    
    }       
      
    
class Stream{

    // general stream model parameters/constraints
    int numStreams, numReachMax;
    double x0, xf, y0, yf, z0;
    double L0, slopeMean, slopeStdev, orientStdev, wMin, wMax, rankAtten;
    double n0, K0, junctionMin;
    
    // constructor
    this(){
        string lineInput;
        string[][] parsedInput; 
        auto inputFile = File("stream.txt","r");                // read in contents of input file
        while (!inputFile.eof()) {
            lineInput = inputFile.readln();
            parsedInput ~= split(lineInput);
            }
        inputFile.close();          
        this.numStreams = to!int(parsedInput[0][1]);            // number of streams/drainages
        this.numReachMax = to!int(parsedInput[1][1]);           // number of reaches; main stream
        this.rankAtten = to!double(parsedInput[2][1]);          // stream length and angle attenuation with rank
        this.x0 = to!double(parsedInput[3][1]);                 // default boundaries (will adjust to match stream footprints)
        this.xf = to!double(parsedInput[4][1]);
        this.y0 = to!double(parsedInput[5][1]);
        this.yf = to!double(parsedInput[6][1]);
        this.z0 = to!double(parsedInput[7][1]);                 // reference elevation at outlet
        this.L0 = to!double(parsedInput[8][1]);                 // length of individual reach sections
        this.slopeMean = to!double(parsedInput[9][1]);          // mean stream bed slope
        this.slopeStdev = to!double(parsedInput[10][1]);        // standard deviation of stream bed slope
        this.orientStdev = to!double(parsedInput[11][1]);       // stream orientation standard deviation between reaches
        this.wMin = to!double(parsedInput[12][1]);              // min & max stream widths (inversely proportional to slope)
        this.wMax = to!double(parsedInput[13][1]); 
        this.n0 = to!double(parsedInput[14][1]);                // default Manning coefficient
        this.K0 = to!double(parsedInput[15][1]);                // default hydraulic conductivity
        this.junctionMin = to!double(parsedInput[16][1]);       // minimum separation distance of stream junctions
        writeln("Read stream model parameters.");       
        }
    
        
    bool CheckCloseStream(double x, double y, double[] xSet, double[] ySet){
        // utility method to check proximity of reach to locations defined by xSet & ySet
        bool closeStream = false;
        for (int i = 1; i < xSet.length; ++i){
            if (Distance(x, y, xSet[i], ySet[i]) <= junctionMin){
                closeStream = true;
                break;
                }
            } 
        return closeStream;
        }

    
    Reach[] RandomWalk(int indexStart, Reach[] reach){
    
        // create set of reaches comprising a stream
        double xStart, yStart, zStart, zEnd;
        double x, y, z, w, orient, slope, orient0;
        int rank, numLocalReach;
        bool fixed, closeStream;
        int numReach = reach.length;        
        
        static double[] xJunction;          // track old stream junction locations to facilitate separation
        static double[] yJunction;
        
        // choose starting location
        if (indexStart == -1){
            // initialize the main stream
            xStart = xf;
            yStart = BoundedRandUniform(yf/2., yf/5., 0., yf);      // pick close to middle
            zStart = z0;
            orient = PI;
            fixed = true;                               // exit point for main stream
            rank = 0;
            }
        else {
            // initialize a tributary stream
            xStart = reach[indexStart].x;
            yStart = reach[indexStart].y;            
            closeStream = CheckCloseStream(xStart, yStart, xJunction, yJunction);
            while (closeStream == true){
                indexStart = uniform(1, reach.length);
                xStart = reach[indexStart].x;
                yStart = reach[indexStart].y;                 
                closeStream = CheckCloseStream(xStart, yStart, xJunction, yJunction);
                }
            zStart = reach[indexStart].z;
            orient0 = reach[indexStart].orient;
            if (orient0 <= PI){orient = orient0 + PI/3.;}
            else {orient = orient0 - PI/3.;}
            fixed = false;
            rank = reach[indexStart].rank + 1;
            }

        // add stream initiation point to list
        xJunction ~= xStart;
        yJunction ~= yStart;  
            
        // first reach of this stream
        x = xStart + 0.5 * (L0*cos(orient));
        y = yStart + 0.5 * (L0*sin(orient));
        slope = BoundedRandUniform(slopeMean, slopeStdev, 0., 1.);
        zEnd = zStart + L0*slope;
        z = 0.5*zStart + 0.5*zEnd;            
        w = wMin + CumNorm(slope, slopeMean, slopeStdev)*(wMax-wMin);
        reach ~= new Reach(x, y, z, L0, w, orient, n0, K0, fixed, rank, zStart, zEnd);
        numReach++;
    
        // random walk to maximum number of reaches for this rank
        numLocalReach = to!int(numReachMax * rankAtten^^rank);
        for (int i = 1; i < numLocalReach; ++i){
            orient = BoundedRandUniform(reach[numReach-1].orient, orientStdev * rankAtten^^rank, PI/2., 3.*PI/2.);
            x = reach[numReach-1].x + 0.5 * (L0*cos(orient));
            y = reach[numReach-1].y + 0.5 * (L0*sin(orient));
            slope = BoundedRandUniform(slopeMean, slopeStdev, 0., 1.);
            zStart = reach[numReach-1].zEnd;
            zEnd = zStart + L0*slope;
            z = 0.5*zStart + 0.5*zEnd;            
            w = wMin + CumNorm(slope, slopeMean, slopeStdev)*(wMax-wMin);            
            fixed = false;
            reach ~= new Reach(x, y, z, L0, w, orient, n0, K0, fixed, rank, zStart, zEnd);
            numReach++;
            }
        
        // process connections and write to connections file
        Connections(reach, numReach, numLocalReach, indexStart);
        
        return reach;
        }

        
    void Connections(Reach[] reach, int numReach, int numLocalReach, int indexStart){
        // process connections between reaches
        int[] up;
        int[] down;
        if (indexStart == -1){
            // main stream
            for (int i = 1; i < numReach; ++i){
                down ~= i-1;
                up ~= i;
                }
            }
        else {
            // tributary
            up ~= numReach-numLocalReach;
            down ~= indexStart;
            for (int i = numReach-numLocalReach+1; i < numReach; ++i){
                down ~= i-1;
                up ~= i;
                }
            }            
        // append to connections file
        auto outputFile = File("connections.csv", "a"); 
        for (int i = 0; i < up.length; ++i){
            string outputLine = to!string(up[i]) ~ "," ~ to!string(down[i]);
            outputFile.writeln(outputLine);
            }
        outputFile.close(); 
        }
        
    }    
    
    
class GridCell: Point{
    // surface grid cell (x-y location, elevation, list of reaches, and flag indicating isolation from stream)
    bool locked;
    int[] reachList;
    this(double x, double y, double z, int[] reachList){
        super(x, y, z);                 // GridCell class inherits location and elevation from Point class
        this.reachList = reachList;
        this.locked = false;
        }
    }
        
        
class Surface{

    // land surface grid and points set, exhibiting erosion associated with stream network
    int intDense, intStream;
    double zStart, dzLevel;
    double gridWt, x0, xf, y0, yf, dx;
    Tuple!(double, double, double, double) bounds;
    Tuple!(double, double) elev;
    GridCell[] gridCell;
    Point[] rPoint;
    
    // constructor
    this(Stream stream, Reach[] reach){
    
        int idCell;
        double zFrac, dzFrac;       
        int[] reachList;
        string lineInput;
        string[][] parsedInput;        
        auto inputFile = File("terrain.txt","r");           // read in contents of input file
        while (!inputFile.eof()) {
            lineInput = inputFile.readln();
            parsedInput ~= split(lineInput);
            }
        inputFile.close();          
        this.intDense = to!int(parsedInput[0][1]);          // grid interval density (must be multiples of 2)
        this.gridWt = to!double(parsedInput[1][1]);         // weighting factor assigned to grid centers (versus stream centers)
        
        // determine starting elevation and elevation increment, per grid level refinement step        
        zFrac = to!double(parsedInput[2][1]);               // starting land elevation, relative to stream network delta z
        dzFrac = to!double(parsedInput[3][1]);              // land elevation uptick per level, relative to stream network delta z       
        elev = SetElevs(zFrac, dzFrac, reach);
        this.zStart = elev[0]; 
        this.dzLevel = elev[1];
        
        bounds = BoundingSquare(stream, reach);
        this.x0 = bounds[0];                                // boundaries of (square) terrain
        this.xf = bounds[1];        
        this.y0 = bounds[2];
        this.yf = bounds[3];        
        this.dx = (xf-x0)/intDense;                         // grid spacing

        // define grid cell objects (as stacked 1-D array) and assign initial elevations
        for (int j = 0; j < intDense; ++j){
            for (int i = 0; i < intDense; ++i){             // fill out by row, then stack rows
                this.gridCell ~= new GridCell(x0+(i+0.5)*dx, y0+(j+0.5)*dx, zStart, reachList);
                }
            }
        this.rPoint = [];       // rPoint = set of points randomly sampled from gridCell
        
        // assign respective arrays of stream reaches to effected cells
        for (int i = 0; i < reach.length; ++i){
            idCell = indexCell(reach[i].x-x0, reach[i].y-y0);
            gridCell[idCell].reachList ~= i;
            }
        
        writeln("Processed terrain model parameters.");       
        }


    Tuple!(double, double) SetElevs(double zFrac, double dzFrac, Reach[] reach){
        // normalize elevation parameters relative to delta z across stream network
        double dzReach, zReachMax;
        double[] zr;
        foreach (rch; reach){zr ~= rch.z;}
        zReachMax = zr.maxElement;                      // maximum elevation in stream network
        dzReach = zr.maxElement - zr.minElement;        // total elevation span of stream network
        zStart = zReachMax + zFrac*dzReach;
        dzLevel = dzFrac * dzReach;
        return tuple(zStart, dzLevel);
        }
        
    
    int indexCell(double x, double y){
        // find index number of grid cell corresponding to location (x, y)
        int iRow, iCol;
        iCol = to!int(x/dx);
        iRow = to!int(y/dx);
        return iRow*intDense + iCol;
        }

        
    int[][] Quad(int level){
        // sets of gridCell indices associated with level no. of binary divisions 
        int iQuadCol, iQuadRow, indexQuad;
        int[][] quadSets;                                   // quadSet[section][list of indices]        
        int numDivide = to!int(2^^level);                   // number of divisions along terrain boundary
        int quadCellSize = intDense/numDivide;              // number of grid cells along each quad section
        quadSets.length = numDivide^^2;
        for (int j = 0; j < intDense; ++j){
            for (int i = 0; i < intDense; ++i){             // by row, then stack rows
                iQuadCol = to!int(i/quadCellSize);          // determine quad section index
                iQuadRow = to!int(j/quadCellSize);
                indexQuad = iQuadRow*numDivide + iQuadCol;
                quadSets[indexQuad] ~= j*intDense + i;      // assign gridCell index (right side) to correct quad section
                }        
            }
        return quadSets;
        }
       
    
    Tuple!(double, double, double, double) BoundingSquare(Stream stream, Reach[] reach){
        // define terrain boundaries; must be square to divide into nested quads
        double x0, xf, y0, yf, xExt, yExt, xm, ym, span;
        double[] xr;
        double[] yr;
        foreach(rch; reach){
            xr ~= rch.x;
            yr ~= rch.y;            
            }
        x0 = min(stream.x0, xr.minElement);           // boundaries possibly extended by extreme positions of stream reach centers
        xf = max(stream.xf, xr.maxElement);
        y0 = min(stream.y0, yr.minElement);
        yf = max(stream.yf, yr.maxElement);            
        xExt = xf - x0;
        yExt = yf - y0;
        xm = (x0+xf)/2.;
        ym = (y0+yf)/2.;        
        span = 1.01*max(xExt, yExt);                 // side of domain, + 1%; scale orthogonal direction to create square
        x0 = xm - span/2.;
        xf = xm + span/2.;
        y0 = ym - span/2.;
        yf = ym + span/2.;        
        return tuple(x0, xf, y0, yf);
        }    

        
    void WtElev(int[] cellList, Reach[] reach, int iLevel){
        // update weighted elevations in subset of grid cells included in cellList
        double xStart, xEnd, yStart, yEnd, xp, yp, zp;
        double totStream = 0.;
        int numStream = 0;
        foreach (indexCell; cellList){
            foreach (indexReach; gridCell[indexCell].reachList){
                totStream += reach[indexReach].z;
                numStream += 1;
                }
            }
        if (numStream > 0){
            foreach (indexCell; cellList){
                gridCell[indexCell].z = (cellList.length*gridWt*gridCell[indexCell].z 
                    + totStream)/(cellList.length*gridWt + numStream);
                }
            }
        else{
            if (gridCell[cellList[0]].locked==false){
                // spatial extent of quad section
                xStart = gridCell[cellList.minElement].x - 0.5*dx;
                xEnd = gridCell[cellList.maxElement].x + 0.5*dx;                
                yStart = gridCell[cellList.minElement].y - 0.5*dx;
                yEnd = gridCell[cellList.maxElement].y + 0.5*dx; 

                // pick random points for subsequent interpolation
                for (int i = 0; i < 2; ++i){
                    xp = xStart + uniform01()*(xEnd-xStart);
                    yp = yStart + uniform01()*(yEnd-yStart);
                    zp = gridCell[cellList[0]].z + dzLevel*(uniform01()/iLevel)^^2.;
                    rPoint ~= new Point(xp, yp, zp);
                    }
                
                // mark cells in this quad as locked so addition points won't be picked subsequently
                foreach (indexCell; cellList){gridCell[indexCell].locked = true;}                
                }
            }
        }    
        
    
    void ShapeTerrain(Reach[] reach){
        // adjust local terrain elevation through progressive quadrant weighted-averaging scheme 
        int level;
        int[][] quadSets;
        level = to!int(log(intDense)/log(1.999));       // number of levels with which to ultimately (binary) subdivide domain
        for (int i = 1; i <= level; ++i){
            writeln("\tProcessing quad level " ~ "\t" ~ to!string(i));
            quadSets = Quad(i);                     // divide domain into quads
            foreach (cellList; quadSets){           // update elevations in each individual quad section and pick points
                WtElev(cellList, reach, i);
                }
            }
        }

        
    void ProcessOutlet(Reach[] reach, Stream stream){
        // (1) remove all random points beyond original stream boundary
        for (int i = rPoint.length-1; i == 0; --i){
            if (rPoint[i].x > stream.xf){remove(rPoint, i);}      // remove point beyond stream system exit
            }
        // (2) add precise stream outlet location to random point set
        rPoint ~= new Point(stream.xf, reach[0].y, stream.z0);
        }

        
    void FillStream(Reach[] reach, Stream stream){
        // add points along stream reaches to random point set, for weighting
        double xrUp, yrUp, drz;
        double intStream = to!int(stream.L0/(2*dx));       // number of intervals along a reach
        double dStream = stream.L0/intStream;              // step size along reach
        foreach(rch; reach){
            drz = (rch.zEnd-rch.zStart)/intStream;
            xrUp = rch.z + 0.5*stream.L0*cos(rch.orient);
            yrUp = rch.z + 0.5*stream.L0*sin(rch.orient);
            for (int i = 0; i < intStream; ++i){                    // add contouring points, stepping back along reach
                rPoint ~= new Point(xrUp + i*dx*cos(-rch.orient),
                yrUp + i*dx*sin(-rch.orient),
                rch.zEnd - i*drz);
                }
            }
        }
 
    
    void WritePoints(Reach[] reach){
        // write random points on surface to output file
        string pointsLine;
        auto outputFile = File("points.csv", "w"); 
        string header = "x" ~ "," ~ "y" ~ "," ~ "z";
        outputFile.writeln(header);
        foreach (pt; rPoint){
            pointsLine = to!string(pt.x) ~ "," ~ to!string(pt.y) ~ "," ~ to!string(pt.z);
            outputFile.writeln(pointsLine);        
            }
        outputFile.close();
        writeln("Wrote random points to file.");    
        }
    
    }
   
        
void main(){

    int indexStart = -1;
    Stream stream;
    Reach[] reach;
    Surface surface;

    // set up stream object
    stream = new Stream();
    
    // write header for connections file
    auto outputFile = File("connections.csv", "w"); 
    string header = "upstream" ~ "," ~ "downstream";
    outputFile.writeln(header); 
    outputFile.close(); 
    
    // step through numStreams and populate reaches
    for (int i = 0; i < stream.numStreams; ++i){
        writeln("\tProcessing stream " ~ "\t" ~ to!string(i));
        reach = stream.RandomWalk(indexStart, reach);
        indexStart = uniform(1, reach.length);
        }
    
    // write reaches file
    WriteReaches(reach);
    
    // process associated land surface
    surface = new Surface(stream, reach);
    surface.ShapeTerrain(reach);                // create guide surface and random sampled points
    surface.FillStream(reach, stream);          // create additional points along streams (for contouring)
    surface.ProcessOutlet(reach, stream);       // truncate points along right boundary at stream outlet    
    surface.WritePoints(reach);                 // write surface to .csv file
    
    writeln("Finished.");
    
    }    
