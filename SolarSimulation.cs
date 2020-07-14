using System;
using System.Collections.Generic;
using System.IO;
using System.Text;


namespace SolarPanelSimulation
{
    //assumptions:
    //all radiation bodies are gray and opaque
    //diameters of pipes are the same throughout the system
    //flow of air around the system should be at least 0.1 so that free convection can still be reasonably modeled
    //since the fluid is water and less heat exchange occurs in first pipe, flows are taken to be fully developed
    //SI units are used in calculations and constants, so construction parameters must be metric
    //all fields should be provided as positive nonzero values
    //emissivity values must be betwen zero and one
    //any temperatures should be given in Kelvin, between 273 and 373 to ensure fluid and thermal properties can be located on tables
    //constants and non-dimensional numbers come from the NCEES FE reference handbook

        /// <summary>
        /// This class models a system in which pumped fluid is exposed to sunlight and ambient air in order to absorb heat from the sunlight while also losing heat to the surroundings.
        /// This model simulates the heating a flat pipe structure that has a nonzero wall thickness throughout by modeling differential sections of the pipe as small pipes with constant heat flux through the surface into the fluid.
        /// By iterating over a multitude of small sections of the pipe, a system that would usually prove difficult to model mathematically without advanced knowledge of fluid and heat mechanics can be simulated albeit with inevitable errors.
        /// @author Parker Beckett
        /// </summary>
    class SolarSimulation
    {
        // fields will be defined in constructor comment
        double pipeToPanel;
        double pipeInPanel;
        double pipeInnerDiameter;
        double pipeOuterDiameter;
        double airTemperature;
        double waterBulkTemperature;
        double kOfPipe;
        double kOfPanel;
        double pipeEmissivity;
        double panelEmissivity;
        double windSpeed;
        double sunScale;
        double solarPower;
        double pumpFlowRate;
        double unitLength;
        bool viewData;

        //constants
        static double sigma = 5.67 * Math.Pow(10, -8);        
        static double CpAir = 1006; //Cp of air is relatively constant up to 50 C, which is close to the record hottest day and thus temperatures above it are unlikely

        //storage of air thermal properties
        static Dictionary<double, double> muTableAir= new Dictionary<double, double>();
        static Dictionary<double, double> nuTableAir= new Dictionary<double, double>();
        static Dictionary<double, double> kTableAir = new Dictionary<double, double>();

        //storage of water thermal properties
        static Dictionary<double, double> muTableH2O = new Dictionary<double, double>();
        static Dictionary<double, double> nuTableH2O = new Dictionary<double, double>();
        static Dictionary<double, double> CpTableH2O = new Dictionary<double, double>();

        /// <summary>
        /// Constructor class that allows the user to define important fields. Also calls the method used to construct property tables.
        /// </summary>
        /// <param name="pipeToPanel">Length of pipe outside of solar panel</param>
        /// <param name="pipeInPanel">Length of pipe inside of solar panel</param>
        /// <param name="pipeInnerDiameter">Inner diameter of pipe</param>
        /// <param name="pipeOuterDiameter">Outer diameter of pipe</param>
        /// <param name="airTemperature">Ambient air temperature</param>
        /// <param name="waterBulkTemperature">Bulk temperature of water in system</param>
        /// <param name="kOfPipe">Thermal conductivity of pipe outside panel</param>
        /// <param name="kOfPanel">Thermal conductivity of pipe inside panel</param>
        /// <param name="pipeEmissivity">Emissivity of pipe</param>
        /// <param name="panelEmissivity">Emissivity of panel</param>
        /// <param name="windSpeed">Non-negative wind speed</param>
        /// <param name="sunScale">Scale from zero to one to simulate weather</param>
        /// <param name="solarPower">Energy output of the sun</param>
        /// <param name="pumpFlowRate">Volume flow rate of system</param>
        /// <param name="unitLength">Differential length to improve accuracy</param>
        /// <param name="viewData">Allows user to view data for each incremental calculation</param>
        SolarSimulation(double pipeToPanel, double pipeInPanel, double pipeInnerDiameter, double pipeOuterDiameter, double airTemperature, double waterBulkTemperature, double kOfPipe, double kOfPanel, double pipeEmissivity, double panelEmissivity, double windSpeed, double sunScale, double solarPower, double pumpFlowRate, double unitLength, bool viewData)
        { 
            this.pipeToPanel = pipeToPanel;
            this.pipeInPanel = pipeInPanel;
            this.pipeInnerDiameter = pipeInnerDiameter;
            this.pipeOuterDiameter = pipeOuterDiameter;
            this.airTemperature = airTemperature;
            this.waterBulkTemperature = waterBulkTemperature;
            this.kOfPipe = kOfPipe;
            this.kOfPanel = kOfPanel;
            this.pipeEmissivity = pipeEmissivity;
            this.panelEmissivity = panelEmissivity;
            this.windSpeed = windSpeed;
            this.sunScale = sunScale;
            this.solarPower = solarPower;
            this.pumpFlowRate = pumpFlowRate;
            this.unitLength = unitLength;
            this.viewData = viewData;
            BuildPropertyTables();
        }

        /// <summary>
        /// Constructs tables of air and water properties for fast and easy access. Accessing separate tables for each property is much faster than searching a  single large table.
        /// </summary>
        public void BuildPropertyTables()
        {
            string lineA = "";
            try
            {
                string dir = Directory.GetCurrentDirectory();
                string path = Path.Combine(dir, "TempMuNuKAir.txt");
                StreamReader r = new StreamReader(path);
                while (lineA != null)
                {
                    lineA = r.ReadLine();
                    string[] elements = lineA.Split(" ");
                    double[] numbers = new double[elements.Length];
                    for (int i = 0; i < elements.Length; i++)
                    {
                        numbers[i] = double.Parse(elements[i]);
                    }
                    muTableAir.Add(numbers[0], numbers[1]);
                    nuTableAir.Add(numbers[0], numbers[2]);
                    kTableAir.Add(numbers[0], numbers[3]);
                }
                r.Close();
            }
            catch (Exception e)
            {
            }
            string lineW = "";
            try
            {
                string dir = Directory.GetCurrentDirectory();
                string path = Path.Combine(dir, "TempMuNuCpH2O.txt");
                StreamReader r = new StreamReader(path);
                while (lineW != null)
                {
                    lineW = r.ReadLine();
                    string[] elements = lineW.Split(" ");
                    double[] numbers = new double[elements.Length];
                    for (int i = 0; i < elements.Length; i++)
                    {
                        numbers[i] = double.Parse(elements[i]);
                    }
                    muTableH2O.Add(numbers[0], numbers[1]);
                    nuTableH2O.Add(numbers[0], numbers[2]);
                    CpTableH2O.Add(numbers[0], numbers[3]);
                }
                r.Close();
            }
            catch (Exception e)
            {
            }
        }

        /// <summary>
        /// Accesses property tables and interpolates accurate values.
        /// </summary>
        /// <param name="property"></param>
        /// <param name="temp"></param>
        /// <returns></returns>
        public double LookupProperty(string substance, string property, double temp)
        {
            Dictionary<double, double> propertyTable = new Dictionary<double, double>();
            if (substance == "air")
            {
                if (property == "mu")
                    propertyTable = muTableAir;
                else if (property == "nu")
                    propertyTable = nuTableAir;
                else if (property == "k")
                    propertyTable = kTableAir;
            }
            else
            {
                if (property == "mu")
                    propertyTable = muTableH2O;
                else if (property == "nu")
                    propertyTable = nuTableH2O;
                else if (property == "Cp")
                    propertyTable = CpTableH2O;
            }

            //interpolates between consecutive table entries
            int x1 = (int)temp / 10 * 10 + 3;
            int x2 = x1 + 10;
            double y1 = propertyTable[x1];
            double y2 = propertyTable[x2];
            return (temp - (double)x1) * ((y2 - y1) / ((double)x2 - (double)x1)) + y1;
        }

        /// <summary>
        /// Calculates Reynolds number based on geometric, and flow properties.
        /// </summary>
        /// <param name="veloc">Flow velocity</param>
        /// <param name="sigDim">Significant dimension</param>
        /// <param name="nu">Kinematic viscosity from tables</param>
        /// <returns>Reynolds number for given parameters</returns>
        public static double ReynoldsNumber(double veloc, double sigDim, double nu)
        {
            return veloc * sigDim / nu;
        }

        /// <summary>
        /// Calculates Prandtl number based on thermal and flow properties.
        /// </summary>
        /// <param name="Cp">Specific heat of substance</param>
        /// <param name="mu">Dynamic viscosity from table</param>
        /// <param name="k">Thermal conductivity from table</param>
        /// <returns>Prandtl number for different parameters</returns>
        public static double PrandtlNumber(double Cp, double mu, double k)
        {
            return Cp * mu / k;
        }

        /// <summary>
        /// Determines the density of water at a particular temperature by using the ratio of its viscosities.
        /// </summary>
        /// <param name="temperature">Temperature of the water</param>
        /// <returns>Density of the water at this temperature</returns>
        public double WaterDensity(double temperature)
        {
            return LookupProperty("H2O", "mu", temperature) / LookupProperty("H2O", "nu", temperature);
        }

        /// <summary>
        /// Solves for the thermal convection coefficient of a forced convection system involving a cylinder. 
        /// This approach is preferable to free convection as the nusselt number requires no additional constants.
        /// </summary>
        /// <param name="Re">Reynolds number of fluid</param>
        /// <param name="Pr">Prandtl number of fluid</param>
        /// <param name="signifDim">Dimension normal to flow, in this case diameter</param>
        /// <param name="k">Thermal conductivity from table</param>
        /// <returns>Convection coefficient of fluid</returns>
        public static double SolveHBarExternal(double Re, double Pr, double signifDim, double k)
        {
            double C;
            double n;

            if (Re < 4)
            {
                C = 0.989;
                n = 0.330;
            }
            else if (Re < 40)
            {
                C = 0.911;
                n = 0.385;
            }
            else if (Re < 4000)
            {
                C = 0.683;
                n = 0.466;
            }
            else if (Re < 40000)
            {
                C = .193;
                n = .618;
            }
            else
            {
                C = 0.0266;
                n = 0.805;
            }

            double Nu = C * Math.Pow(Re, n) * Math.Pow(Pr, 0.33);
            return Nu * k / signifDim;
        }

        /// <summary>
        /// Uses the Newton-Raphson method for finding higher-order roots to ensure accurate values. Coefficients have been obtained by hand, with values being substituted here to save computing power.
        /// </summary>
        /// <param name="Ts">Surface temperature of cylinder</param>
        /// <param name="emissivity">Emissivity of the material</param>
        /// <param name="hBarExt">Convective coefficient for flow around system</param>
        /// <param name="kSubst">Conductivity coefficient for material</param>
        /// <returns></returns>
        public double NewtonRaphsonSolver(double Ts, double emissivity, double hBarExt, double kSubst)
        {
            double perimeter = Math.PI * pipeOuterDiameter;
            double coeffA = (1-emissivity) * solarPower * sunScale * perimeter * unitLength / 2; 
            double coeffB = hBarExt * perimeter * unitLength * airTemperature;
            double coeffC = (2 * Math.PI * kSubst * unitLength * waterBulkTemperature) / Math.Log(pipeOuterDiameter/ pipeInnerDiameter);
            double coeffD = emissivity * sigma * perimeter * unitLength / 2;            
            double coeffE = coeffB / airTemperature;
            double coeffF = coeffC / waterBulkTemperature;
            double correction = SolveTs(Ts, coeffA, coeffB, coeffC, coeffD, coeffE, coeffF) / SolveTsDeriv(Ts, coeffD, coeffE, coeffF);
            while (Math.Abs(correction) > 0.00001)
            {
                correction = SolveTs(Ts, coeffA, coeffB, coeffC, coeffD, coeffE, coeffF) / SolveTsDeriv(Ts, coeffD, coeffE, coeffF);
                Ts -= correction;
            }
            return Ts;
        }

        /// <summary>
        /// Helper method that compacts equation balancing, defines actual function.
        /// </summary>
        /// <param name="Ts">Surface temperature of cylinder</param>
        /// <param name="coeffA"></param>
        /// <param name="coeffB"></param>
        /// <param name="coeffC"></param>
        /// <param name="coeffD"></param>
        /// <param name="coeffE"></param>
        /// <param name="coeffF"></param>
        /// <returns></returns>
        public double SolveTs(double Ts, double coeffA, double coeffB, double coeffC, double coeffD, double coeffE, double coeffF)
        {
            return coeffA + coeffB + coeffC - Ts * ((coeffD * Ts * Ts * Ts) + coeffE + coeffF);
        }

        /// <summary>
        /// Helper method that compacts equation balancing, defines derivative function.
        /// </summary>
        /// <param name="Ts">Surface temperature of cylinder</param>
        /// <param name="coeffA"></param>
        /// <param name="coeffB"></param>
        /// <param name="coeffC"></param>
        /// <param name="coeffD"></param>
        /// <param name="coeffE"></param>
        /// <param name="coeffF"></param>
        /// <returns></returns>
        public double SolveTsDeriv(double Ts, double coeffD, double coeffE, double coeffF)
        {
            return -1 * (4 * coeffD * Ts * Ts * Ts * +coeffE + coeffF);
        }

        /// <summary>
        /// Simulates an unsteady-state system by using differential lengths of pipe to analyze thermal behavior.
        /// </summary>
        public void SimulateThermodynamics()
        {
            double flowVelocity = pumpFlowRate / (Math.PI * Math.Pow((pipeInnerDiameter / 2), 2));
            double airReynolds = ReynoldsNumber(windSpeed, pipeOuterDiameter, LookupProperty("air", "nu", airTemperature));
            double airPrandtl = PrandtlNumber(CpAir, LookupProperty("air", "mu", airTemperature), LookupProperty("air", "k", airTemperature));
            double hBarAir = SolveHBarExternal(airReynolds, airPrandtl, pipeOuterDiameter, LookupProperty("air", "k", airTemperature));
            double pos = 0;
            double maxTemp = 0;
            double interTemp;
            bool maxReached = false;

            if(viewData)
            Console.WriteLine("\n **Initial Path** \n");
            for (double i = 0; i < pipeToPanel; i += unitLength)
            {
                interTemp = DetermineUnitTempChange(flowVelocity, unitLength, hBarAir, kOfPipe, pipeEmissivity, maxTemp);
                if (interTemp > maxTemp)
                    maxTemp = interTemp;
                pos += unitLength;
            }
            pos = 0;

            if (viewData)
            Console.WriteLine("\n **Panel Path** \n");
            for (double i = 0; i < pipeInPanel; i += unitLength)
            { 
                interTemp = DetermineUnitTempChange(flowVelocity, unitLength, hBarAir, kOfPanel, panelEmissivity, maxTemp);
                if (interTemp > maxTemp)
                    maxTemp = interTemp;
                else if (interTemp == maxTemp)
                    maxReached = true;
                pos += unitLength;
            }
            pos = 0;

            if (viewData)
                Console.WriteLine("\n **Return Path** \n");
            for (double i = 0; i < pipeToPanel; i += unitLength)
            {
                interTemp = DetermineUnitTempChange(flowVelocity, unitLength, hBarAir, kOfPipe, pipeEmissivity, maxTemp);
                if (interTemp > maxTemp)
                    maxTemp = interTemp;
                pos += unitLength;
            }
            Console.WriteLine();
            Console.WriteLine($"Maximum temperature reached under these conditions? {maxReached}");
            Console.WriteLine($"Final Temperature: {waterBulkTemperature}");
            Console.WriteLine($"Maximum Temperature: {maxTemp}");
        }

        /// <summary>
        /// Determines incremental temperature changes for each differential pipe section.
        /// </summary>
        /// <param name="flowVelocity">Fluid velocity in system</param>
        /// <param name="hBarAir">Convection coefficient of air</param>
        /// <param name="kSubst">Conductivity coefficient of material</param>
        /// <param name="emissivity">Emissivity of material</param>
        public double DetermineUnitTempChange(double flowVelocity, double differentialLength, double hBarAir, double kSubst, double emissivity, double maxTemp)
        {
            double surfaceTemperature = NewtonRaphsonSolver(waterBulkTemperature, emissivity, hBarAir, kSubst);
            double qCond = SolveConductionPerUnitLength(surfaceTemperature, kSubst);
            waterBulkTemperature = SimulateConstantFluxFlow(qCond);
            
            if(viewData)
                Console.WriteLine($"Surface Temperature: {Math.Round(surfaceTemperature, 5)}, Conduction Heat: {Math.Round(qCond, 5)}, Water Temperature: {Math.Round(waterBulkTemperature, 5)}");
            return waterBulkTemperature;
        }

        /// <summary>
        /// Solves for conductive heat transfer based the temperature difference of the previous iteration.
        /// </summary>
        /// <param name="surfaceTemp">Temperature of pipe</param>
        /// <param name="waterTemp">Temperature of water in pipe</param>
        /// <returns></returns>
        public double SolveConductionPerUnitLength(double surfaceTemp, double kSubst)
        {
            return (2 * Math.PI * kSubst * unitLength * (surfaceTemp - waterBulkTemperature)) / Math.Log(pipeOuterDiameter / pipeInnerDiameter);
        }

        /// <summary>
        /// Models each unit section of pipe as a tube with constant surface flux so that the model can be simplified.
        /// </summary>
        /// <param name="qCond">Conduction heat through the tube inner surface</param>
        /// <returns>Change in temperature due to the conductive flux</returns>
        public double SimulateConstantFluxFlow(double qCond)
        {
            double perimeter = 2 * Math.PI * pipeInnerDiameter;
            double conductiveFlux = qCond / (perimeter* unitLength);
            return (conductiveFlux * perimeter) / (pumpFlowRate * WaterDensity(waterBulkTemperature) * LookupProperty("H2O", "Cp", waterBulkTemperature)) * unitLength + waterBulkTemperature;
        }

        // <summary>
        // Driver function for simulation class. Parameters are defined in this method, either by the user or by default.
        // </summary>
        // <param name = "args" ></ param >
        static void Main(string[] args)
        {
            Console.OutputEncoding = Encoding.UTF8;
            bool viewData = false;
            SolarSimulation simulate;

            Console.WriteLine("Would you like to view data about incremental steps? Enter y or n");
            string ynD = Console.ReadLine();
            if (ynD == "y")
                viewData = true;
            else
                viewData = false;

            Console.WriteLine("Would you like to enter your own parameters? Enter y or n");
            string ynP = Console.ReadLine();

            if (ynP == "y")
            {
 
                Console.WriteLine("Length of pipe from tank to panel (m):");
                Console.WriteLine("Default value: 2");
                double a = double.Parse(Console.ReadLine());

                Console.WriteLine("Length of pipe inside panel (m):");
                Console.WriteLine("Default value: 40");
                double b = double.Parse(Console.ReadLine());

                Console.WriteLine("Inner diameter of pipe (m):");
                Console.WriteLine("Default value: 0.019");
                double c = double.Parse(Console.ReadLine());

                Console.WriteLine("Outer diameter of pipe (m):");
                Console.WriteLine("Default value: 0.02");
                double d = double.Parse(Console.ReadLine());

                Console.WriteLine("Air temperature (K):");
                Console.WriteLine("Default value: 298");
                double e = double.Parse(Console.ReadLine());

                Console.WriteLine("Initial water temperature (K):");
                Console.WriteLine("Default value: 288");
                double f = double.Parse(Console.ReadLine());

                Console.WriteLine("Thermal conductivity of pipe to panel:");
                Console.WriteLine("Default value: 100");
                double g = double.Parse(Console.ReadLine());

                Console.WriteLine("Thermal conductivity of pipe in panel:");
                Console.WriteLine("Default value: 400");
                double h = double.Parse(Console.ReadLine());

                Console.WriteLine("Emissivity of pipe to panel:");
                Console.WriteLine("Default value: 0.5");
                double i = double.Parse(Console.ReadLine());

                Console.WriteLine("Emissivity of pipe in panel:");
                Console.WriteLine("Default value: 0.1");
                double j = double.Parse(Console.ReadLine());

                Console.WriteLine("Speed of wind surrounding system (m/s):");
                Console.WriteLine("Default value: 1.0");
                double k = double.Parse(Console.ReadLine());

                Console.WriteLine("Ratio of sun reaching panel (0-1):");
                Console.WriteLine("Default value: 0.75");
                double l = double.Parse(Console.ReadLine());

                Console.WriteLine("Solar intensity (W/m^2):");
                Console.WriteLine("Default value: 1353");
                double m = double.Parse(Console.ReadLine());

                Console.WriteLine("Flow rate of water in system (m^3/s):");
                Console.WriteLine("Default value: 0.0000033");
                double n = double.Parse(Console.ReadLine());

                Console.WriteLine("Length of differential tube length (m):");
                Console.WriteLine("Default value: 0.005");
                double o = double.Parse(Console.ReadLine());

                Console.WriteLine("Press Enter to start");
                Console.ReadLine();

                simulate = new SolarSimulation(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, viewData);
            }
            else
            {
                Console.WriteLine("Length of pipe from tank to panel (m): 2");
                Console.WriteLine("Length of pipe inside panel (m): 40");
                Console.WriteLine("Inner diameter of pipe (m): 0.019");
                Console.WriteLine("Outer diameter of pipe (m): 0.02");
                Console.WriteLine("Air temperature (K): 298");
                Console.WriteLine("Initial water temperature (K): 288");
                Console.WriteLine("Thermal conductivity of pipe to panel: 100");
                Console.WriteLine("Thermal conductivity of pipe in panel: 400");
                Console.WriteLine("Emissivity of pipe to panel: 0.5");
                Console.WriteLine("Emissivity of pipe in panel: 0.1");
                Console.WriteLine("Speed of wind surrounding system (m/s): 1.0");
                Console.WriteLine("Ratio of sun penetrating clouds (0-1): 0.75");
                Console.WriteLine("Solar intensity (W/m^2): 1353");
                Console.WriteLine("Flow rate of water in system (m^3/s): 0.0000033");
                Console.WriteLine("Length of differential tube length (m): 0.005");

                Console.WriteLine("Press Enter to start");
                Console.ReadLine();

                simulate = new SolarSimulation(2, 40, 0.019, 0.02, 298, 288, 100, 400, 0.5, 0.1, 1.0, 0.75, 1353, 0.0000033, 0.005, viewData);
            }
            simulate.SimulateThermodynamics();

            Console.WriteLine("Would you like to try another simulation? Enter y or n");
            string tryAgain = Console.ReadLine();
            if (tryAgain == "y")
            {
                Console.WriteLine("\n **New Simulation** \n");
                Main(new string[0]);
            }

            Console.WriteLine("\n" + "Thank you for using this software. Consult the ProjectCode.txt file to get a closer look at how the program works.");
            Console.WriteLine("\n" + "\u00a9" + " Parker Beckett 07/2020. Stay safe!" + "\n");
            Console.WriteLine("Press Enter to close");
            Console.ReadLine();

        }
    }
}