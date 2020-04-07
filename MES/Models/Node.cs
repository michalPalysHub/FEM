namespace MES.Models
{
    public class Node
    {
        public int ID;
        public double x, y;  // współrzędne globalne
        public double t;     // temperatura
        public bool bc;      // warunki brzegowe

        public Node(int ID)
        {
            this.ID = ID;
        }

        public Node(int ID, double x, double y, double t)
        {
            this.ID = ID;
            this.x = x;
            this.y = y;
            this.t = t;
        }
    }
}
