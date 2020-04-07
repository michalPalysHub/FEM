using System;
using MathNet.Numerics.LinearAlgebra;

namespace MES.Models
{
    public class ElementUniversal
    {
        // pochodna dN_1/dξ, dN_2/dξ, ... dla każdego punktu całkowania
        public Matrix<double> dN_i_dKsi;

        // pochodna dN_1/dη, dN_2/dη, ... dla każdego punktu całkowania
        public Matrix<double> dN_i_dEta;

        // Wartość N_1, N_2, ... dla każdego punktu całkowania
        public Matrix<double> N_i;


        public int NumberOfNodes = 4;
        public GlobalData globalData;

        public ElementUniversal(int nPc)
        {
            globalData = new GlobalData(nPc);

            dN_i_dKsi = Matrix<double>.Build.Dense(nPc * nPc, NumberOfNodes);
            dN_i_dEta = Matrix<double>.Build.Dense(nPc * nPc, NumberOfNodes);
            N_i = Matrix<double>.Build.Dense(nPc * nPc, NumberOfNodes);

            for (int Pc_i = 0; Pc_i < nPc * nPc; Pc_i++)
            {
                dN_i_dKsi[Pc_i, 0] = -0.25 * (1 - globalData.Pc_2D[Pc_i % nPc]);
                dN_i_dKsi[Pc_i, 1] = 0.25 * (1 - globalData.Pc_2D[Pc_i % nPc]);
                dN_i_dKsi[Pc_i, 2] = 0.25 * (1 + globalData.Pc_2D[Pc_i % nPc]);
                dN_i_dKsi[Pc_i, 3] = -0.25 * (1 + globalData.Pc_2D[Pc_i % nPc]);

                dN_i_dEta[Pc_i, 0] = -0.25 * (1 - globalData.Pc_2D[Pc_i / nPc]);
                dN_i_dEta[Pc_i, 1] = -0.25 * (1 + globalData.Pc_2D[Pc_i / nPc]);
                dN_i_dEta[Pc_i, 2] = 0.25 * (1 + globalData.Pc_2D[Pc_i / nPc]);
                dN_i_dEta[Pc_i, 3] = 0.25 * (1 - globalData.Pc_2D[Pc_i / nPc]);

                N_i[Pc_i, 0] = 0.25 * (1 - globalData.Pc_2D[Pc_i / nPc]) * (1 - globalData.Pc_2D[Pc_i % nPc]);
                N_i[Pc_i, 1] = 0.25 * (1 + globalData.Pc_2D[Pc_i / nPc]) * (1 - globalData.Pc_2D[Pc_i % nPc]);
                N_i[Pc_i, 2] = 0.25 * (1 + globalData.Pc_2D[Pc_i / nPc]) * (1 + globalData.Pc_2D[Pc_i % nPc]);
                N_i[Pc_i, 3] = 0.25 * (1 - globalData.Pc_2D[Pc_i / nPc]) * (1 + globalData.Pc_2D[Pc_i % nPc]);
            }
        }

        public void Print()
        {
            Console.WriteLine("\n\n///////////////////////////////////////////////////");
            Console.WriteLine("////////////////Element uniwersalny////////////////");
            Console.WriteLine("///////////////////////////////////////////////////");
            Console.WriteLine("\nLiczba punktów całkowania: {0}\n", globalData.nPc);

            Console.WriteLine("dN_i/dKsi:");
            Console.WriteLine(dN_i_dKsi.ToString()); 

            Console.WriteLine("dN_i/dEta:");
            Console.WriteLine(dN_i_dEta.ToString());

            Console.WriteLine("N_i:");
            Console.WriteLine(N_i.ToString());
        }
    }
}
