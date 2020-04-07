using System;
using System.Collections.Generic;

namespace MES.Models
{
    public class Element
    {
        public int ID;
        public List<Node> Nodes;
        public List<bool> SidesWithBc;
        public ElementUniversal elementUniversal;

        public Element(int ID, List<Node> Nodes)
        {
            this.ID = ID;
            this.Nodes = Nodes;
            SidesWithBc = new List<bool>(4)
            {
                false, false, false, false
            };
        }

        public void Print()
        {
            Console.WriteLine("ID Elementu = {0}", ID);
            for (int i = 0; i < Nodes.Count; i++)
            {
                Console.WriteLine("\t{0}) ID Węzła = {1}, x = {2}, y = {3}, bc = {4}, t = {5})", i + 1, Nodes[i].ID, Nodes[i].x, Nodes[i].y, Nodes[i].bc, Nodes[i].t);
            }
            Console.WriteLine("bc na bokach:");
            for (int j = 0; j < 4; j++)
            {
                Console.WriteLine("\tbok " + j + ": " + SidesWithBc[j]);
            }
            Console.WriteLine();
        }
    }
}
