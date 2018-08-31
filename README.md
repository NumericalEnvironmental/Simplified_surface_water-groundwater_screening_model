# Simplified_surface_water-groundwater_screening_model

![Preview](https://numericalenvironmental.files.wordpress.com/2018/08/figure_1.jpg)

This python script models the transient flow of water through a network of connected stream reaches of differing slopes and widths in accordance with the Manning equation, also accounting for discharge through the stream bed to underlying groundwater. The solution strategy based upon solving the transient streamflow equations by an explicit finite-difference (forward-differencing in time) scheme. An operator-splitting approach is then implemented to model percolation through the stream bed to underlying groundwater. A second model component uses the streambed discharges as source terms for a groundwater mounding equation (using an analytical expression proposed by Hantush, 1967). The combination of the two models (numerical surface, analytical groundwater) can then be used as a screening tool for engineering calculations or as part of the initial setup for a more complex numerical surface water-groundwater model. Potential applications include managed aquifer recharge by releases of water into otherwise dry stream beds. See my blog post (https://numericalenvironmental.wordpress.com/2018/08/31/a-mixed-numerical-analytical-screening-model-for-groundwater-recharge-from-a-losing-stream/) for more details.

This particular script is really just in alpha test mode at the moment. I have not yet had a chance to compare it to numerical results for a coupled system (e.g., MODFLOW with the SWR package) although the surface water module produces results that are consistent with HEC-RAS, at least for a simple test example. As such, user beware!

This python 3 script requires the numpy, scipy, pandas, and matplotlib libraries. Required text input files include the following:

* aquifer.txt - uniform aquifer hydraulic conductivity, specific storage, thickness, initial head, and vadose zone thickness
* params.txt - various finite difference model constraints; see comments under Params class constructor
* print_times.txt - list of output times
* reaches.csv - reach properties; see comments under Reach class constructor
* sources.txt - reach index, volumetric flow rate, and cessation time for each time-dependent source term
* well.txt - well locations and monitor period (start, end, and number of uniformly-distributed sample points)

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.
THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
