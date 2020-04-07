using System;
using System.Collections.Generic;

namespace MES.Models
{
    public class GlobalData
    {
        // liczba punktów całkowania
        public int nPc;

        // punkty całkowania
        public List<double> Pc_2D;

        public List<double> Pc_1D_Ksi;
        public List<double> Pc_1D_Eta;

        // wagi całkowania
        public List<double> Wc_2D;

        public GlobalData(int nPc)
        {
            this.nPc = nPc;
            Pc_1D_Ksi = new List<double>(4 * nPc);
            Pc_1D_Eta = new List<double>(4 * nPc);
            Pc_2D = new List<double>(nPc);
            Wc_2D = new List<double>(nPc);

            if (nPc == 2)
            {
                Pc_2D.Add((-1) / Math.Sqrt(3));
                Pc_2D.Add(1 / Math.Sqrt(3));

                Wc_2D.Add(1);
                Wc_2D.Add(1);

                // można tutaj zainicjować 2 razy mniej tych punktów i w programie dodawać znak - powyżej 4tego elementu
                Pc_1D_Ksi.Add((-1) / Math.Sqrt(3));
                Pc_1D_Ksi.Add(1 / Math.Sqrt(3));
                Pc_1D_Ksi.Add(1);
                Pc_1D_Ksi.Add(1);
                Pc_1D_Ksi.Add(1 / Math.Sqrt(3));
                Pc_1D_Ksi.Add((-1) / Math.Sqrt(3));
                Pc_1D_Ksi.Add(-1);
                Pc_1D_Ksi.Add(-1);

                Pc_1D_Eta.Add(-1);
                Pc_1D_Eta.Add(-1);
                Pc_1D_Eta.Add((-1) / Math.Sqrt(3));
                Pc_1D_Eta.Add(1 / Math.Sqrt(3));
                Pc_1D_Eta.Add(1);
                Pc_1D_Eta.Add(1);
                Pc_1D_Eta.Add(1 / Math.Sqrt(3));
                Pc_1D_Eta.Add((-1) / Math.Sqrt(3));
            }
            else if (nPc == 3)
            {
                Pc_2D.Add(-0.77);
                Pc_2D.Add(0);
                Pc_2D.Add(0.77);

                Wc_2D.Add(0.5555555555555556);
                Wc_2D.Add(0.8888888888888889);
                Wc_2D.Add(0.5555556666655556);

                Pc_1D_Ksi.Add(-0.77);
                Pc_1D_Ksi.Add(0);
                Pc_1D_Ksi.Add(0.77);
                Pc_1D_Ksi.Add(1);
                Pc_1D_Ksi.Add(1);
                Pc_1D_Ksi.Add(1);
                Pc_1D_Ksi.Add(0.77);
                Pc_1D_Ksi.Add(0);
                Pc_1D_Ksi.Add(-0.77);
                Pc_1D_Ksi.Add(-1);
                Pc_1D_Ksi.Add(-1);
                Pc_1D_Ksi.Add(-1);

                Pc_1D_Eta.Add(-1);
                Pc_1D_Eta.Add(-1);
                Pc_1D_Eta.Add(-1);
                Pc_1D_Eta.Add(-0.77);
                Pc_1D_Eta.Add(0);
                Pc_1D_Eta.Add(0.77);
                Pc_1D_Eta.Add(1);
                Pc_1D_Eta.Add(1);
                Pc_1D_Eta.Add(1);
                Pc_1D_Eta.Add(0.77);
                Pc_1D_Eta.Add(0);
                Pc_1D_Eta.Add(-0.77);
            }
            else
            {
                throw new ArgumentException("Musisz wybrać wariant 2 albo 3 punktów całkowania");
            }
        }
    }
}
