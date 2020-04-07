using System;
using System.Collections.Generic;

namespace MES.Models
{
    public class Grid
    {
        // lista węzłów
        public List<Node> Nodes;

        // lista elementów 
        public List<Element> Elements;

        // wysokość i szerokość siatki
        public double H, W;

        // ilość węzłów na wysokość i na szerokość
        public int nH, nW;

        // ilość węzłów, liczba elementów
        public int nN, nE;


        ////////////////////////////////////////////////
        // zakładamy, że elementy siatki są prostokątne
        public Grid(double H, double W, int nH, int nW, int nPc, double t_0)
        {
            if ((H <= 0) || (W <= 0))
            {
                throw new ArgumentException("Wysokość i szerokość siatki musi byc liczbą dodatnią!");
            }

            if ((nH < 2) || (nW < 2))
            {
                throw new ArgumentException("Siatka musi mieć wymiary co najmniej 2x2!");
            }

            this.H = H;
            this.W = W;
            this.nH = nH;
            this.nW = nW;
            nN = nH * nW;
            nE = (nH - 1) * (nW - 1);

            Nodes = new List<Node>(nN);
            Elements = new List<Element>(nE);


            // tworzenie właściwej siatki
            int column = 0;
            double dH = H / (nH - 1);
            double dW = W / (nW - 1);

            for (int i = 1; i <= nE; i++)
            {
                var n1 = new Node(i + column, column * dW, (i - 1 + column) % nH * dH, t_0);
                var n2 = new Node(i + nH + column, (column + 1) * dW, (i - 1 + column) % nH * dH, t_0);
                var n3 = new Node(i + nH + 1 + column, (column + 1) * dW, (i + column) % nH * dH, t_0);
                var n4 = new Node(i + 1 + column, column * dW, (i + column) % nH * dH, t_0);

                var nodes = new List<Node>
                {
                    n1,
                    n2,
                    n3,
                    n4
                };

                foreach (var node in nodes)
                {
                    if (node.x == 0 || node.x == ((nW - 1) * dW) || node.y == 0 || node.y == ((nH - 1) * dH))
                    {
                        node.bc = true;
                    }
                }

                var element = new Element(i, nodes);
                int numberOfSidesInElement = 4;

                // sprawdzanie, na których ścianach elementu mamy warunki brzegowe
                for (int sideNumber = 0; sideNumber < numberOfSidesInElement; sideNumber++)
                {
                    if (nodes[sideNumber].bc == true && nodes[(sideNumber + 1) % numberOfSidesInElement].bc == true)
                    {
                        element.SidesWithBc[sideNumber] = true;
                    }
                }

                Elements.Add(element);


                // dodawanie węzłów do listy Nodes
                if (Nodes.Count != 0)
                {
                    int numberOfNodesInElement = 4;
                    for (int j = 0; j < numberOfNodesInElement; j++)
                    {
                        bool isInNodes = false;
                        for (int k = 0; k < Nodes.Count; k++)
                        {
                            if (nodes[j].ID == Nodes[k].ID)
                            {
                                isInNodes = true;
                            }
                        }

                        if (isInNodes == false)
                        {
                            Nodes.Add(nodes[j]);
                        }
                    }
                }
                else
                {
                    foreach (var node in nodes) Nodes.Add(node);
                }

                if (i % (nH - 1) == 0) column++;
            }
        }

        public void Print()
        {
            Console.WriteLine("///////////////////////////////////////////////////");
            Console.WriteLine("///////////////////Siatka MES//////////////////////");
            Console.WriteLine("///////////////////////////////////////////////////\n");
            foreach (var Element in Elements)
            {
                Element.Print();
            }
        }
    }
}
