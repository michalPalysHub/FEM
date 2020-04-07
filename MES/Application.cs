using MathNet.Numerics.LinearAlgebra;
using MES.Models;
using System;
using System.Collections.Generic;
using System.IO;

namespace MES
{
    public class Application
    {
        static void Main(string[] args)
        {
            try
            {
                ////////////////////////////////////
                //// Dane globalne dla programu ////
                ////////////////////////////////////

                double H, W, k, alfa, c, ro, t_alfa, t_0;
                int nH, nW, tau, dTau, numberOfIntegrationPoints, numberOfIterationsToPrint;

                using(var reader = File.OpenText("../../data.txt"))
                {
                    H = double.Parse(reader.ReadLine()); // wysokość siatki w metrach
                    W = double.Parse(reader.ReadLine()); // szerokość siatki w metrach
                    nH = int.Parse(reader.ReadLine());     // ilość węzłów na wysokość
                    nW = int.Parse(reader.ReadLine());     // ilość węzłów na szerokość
                    k = double.Parse(reader.ReadLine());     // wartość współczynnika przewodzenia ciepła
                    alfa = double.Parse(reader.ReadLine()); //wartość współczynnika przejmowania ciepła
                    c = double.Parse(reader.ReadLine());    // wartość ciepła właściwego danego materiału
                    ro = double.Parse(reader.ReadLine());  // wartość gęstości danego materiału
                    t_alfa = double.Parse(reader.ReadLine()); // temperatura otoczenia
                    t_0 = double.Parse(reader.ReadLine());     // temperatura początkowa materiału
                    tau = int.Parse(reader.ReadLine());     // czas symulacji
                    dTau = int.Parse(reader.ReadLine());     // krok czasowy
                    
                    numberOfIntegrationPoints = int.Parse(reader.ReadLine());  // liczba punktów całkowania równa 2 lub 3
                    numberOfIterationsToPrint = int.Parse(reader.ReadLine());
                }

                int nIterations = tau / dTau;  // ilość kroków czasowych 

                /*
                double H = 0.1; // wysokość siatki w metrach
                double W = 0.1; // szerokość siatki w metrach
                int nH = 4;     // ilość węzłów na wysokość
                int nW = 4;     // ilość węzłów na szerokość
                double k = 25;     // wartość współczynnika przewodzenia ciepła
                double alfa = 300; //wartość współczynnika przejmowania ciepła
                double c = 700;    // wartość ciepła właściwego danego materiału
                double ro = 7800;  // wartość gęstości danego materiału
                double t_alfa = 1200; // temperatura otoczenia
                double t_0 = 100;     // temperatura początkowa materiału
                int tau = 500;     // czas symulacji
                int dTau = 50;     // krok czasowy
                int nIterations = tau / dTau;  // ilość kroków czasowych 


                int numberOfIntegrationPoints = 3;  // liczba punktów całkowania równa 2 lub 3
                int numberOfIterationsToPrint = 10;
                */

                ///////////////////////////////////////////////////////////////////////////////

                ////////////////////////////////////
                ////// Tworzenie siatki MES ////////
                ////////////////////////////////////

                Grid grid = new Grid(H, W, nH, nW, numberOfIntegrationPoints, t_0);
                //grid.Print();


                ///////////////////////////////////////////////////
                ///// Obliczanie dla każdego elementu siatki: /////
                /////  - macierzy H_l,                        /////
                /////  - macierzy H_bc_l z war.brzeg.,        /////
                /////  - macierzy C_l,                        /////
                /////  - wektora P_l.                         /////
                ///////////////////////////////////////////////////

                // element uniwersalny, z którego będziemy pobierać odpowiednie współrzędne punktów całkowania i wagi całkowania
                var elementUniversal = new ElementUniversal(numberOfIntegrationPoints);
                elementUniversal.Print();
                int numberOfNodesInElement = 4;

                // lista macierzy H_l, H_bc_l, C_l dla każdego elementu siatki
                var H_l_matrices = new List<Matrix<double>>(grid.nE);
                var H_bc_l_matrices = new List<Matrix<double>>(grid.nE);
                var C_l_matrices = new List<Matrix<double>>(grid.nE);
                var P_l_vectors = new List<Matrix<double>>(grid.nE);

                for (int gridElementIndex = 0; gridElementIndex < grid.nE; gridElementIndex++)
                {
                    H_l_matrices.Add(Matrix<double>.Build.Dense(numberOfNodesInElement, numberOfNodesInElement));
                    H_bc_l_matrices.Add(Matrix<double>.Build.Dense(numberOfNodesInElement, numberOfNodesInElement));
                    C_l_matrices.Add(Matrix<double>.Build.Dense(numberOfNodesInElement, numberOfNodesInElement));
                    P_l_vectors.Add(Matrix<double>.Build.Dense(numberOfNodesInElement, 1));
                }

                // dla każdego elementu siatki tworzę macierze H lokalną i C lokalną
                for (int gridElementIndex = 0; gridElementIndex < grid.nE; gridElementIndex++)
                {
                    Element element = grid.Elements[gridElementIndex];

                    // listy macierzy H_i oraz J_i dla każdego punktu całkowania elementu
                    var H_i_matrices = new List<Matrix<double>>(numberOfIntegrationPoints * numberOfIntegrationPoints);
                    var C_i_matrices = new List<Matrix<double>>(numberOfIntegrationPoints * numberOfIntegrationPoints);
                    var J_i_matrices = new List<Matrix<double>>(numberOfIntegrationPoints * numberOfIntegrationPoints);


                    // obliczanie wszystkich H_i, C_i oraz J_i dla danego elementu
                    for (int integrationPointIndex = 0; integrationPointIndex < numberOfIntegrationPoints * numberOfIntegrationPoints; integrationPointIndex++)
                    {
                        var dN_i_dx = Matrix<double>.Build.Dense(numberOfNodesInElement, 1);
                        var dN_i_dy = Matrix<double>.Build.Dense(numberOfNodesInElement, 1);
                        var N_i = Matrix<double>.Build.Dense(numberOfNodesInElement, 1);

                        var J_i_matrix = Matrix<double>.Build.Dense(2, 2);

                        for (int shapeFunctionIndex = 0; shapeFunctionIndex < numberOfNodesInElement; shapeFunctionIndex++)
                        {
                            J_i_matrix[0, 0] += elementUniversal.dN_i_dKsi[integrationPointIndex, shapeFunctionIndex] * element.Nodes[shapeFunctionIndex].x;
                            J_i_matrix[1, 0] += elementUniversal.dN_i_dEta[integrationPointIndex, shapeFunctionIndex] * element.Nodes[shapeFunctionIndex].x;
                            J_i_matrix[0, 1] += elementUniversal.dN_i_dKsi[integrationPointIndex, shapeFunctionIndex] * element.Nodes[shapeFunctionIndex].y;
                            J_i_matrix[1, 1] += elementUniversal.dN_i_dEta[integrationPointIndex, shapeFunctionIndex] * element.Nodes[shapeFunctionIndex].y;
                        }

                        var J_i_matrix_inversed = J_i_matrix.Inverse();

                        for (int j = 0; j < numberOfNodesInElement; j++)
                        {
                            dN_i_dx[j, 0] = J_i_matrix_inversed[0, 0] * elementUniversal.dN_i_dKsi[integrationPointIndex, j] + J_i_matrix_inversed[0, 1] * elementUniversal.dN_i_dEta[integrationPointIndex, j];
                            dN_i_dy[j, 0] = J_i_matrix_inversed[1, 0] * elementUniversal.dN_i_dKsi[integrationPointIndex, j] + J_i_matrix_inversed[1, 1] * elementUniversal.dN_i_dEta[integrationPointIndex, j];
                            N_i[j, 0] = elementUniversal.N_i[integrationPointIndex, j];
                        }

                        Matrix<double> H_i_matrix = k * ((dN_i_dx * dN_i_dx.Transpose()) + (dN_i_dy * dN_i_dy.Transpose()));
                        Matrix<double> C_i_matrix = c * ro * (N_i * N_i.Transpose());

                        H_i_matrices.Add(H_i_matrix);
                        C_i_matrices.Add(C_i_matrix);
                        J_i_matrices.Add(J_i_matrix);
                    }

                    for (int integrationPointIndex = 0; integrationPointIndex < H_i_matrices.Capacity; integrationPointIndex++)
                    {
                        H_l_matrices[gridElementIndex] += H_i_matrices[integrationPointIndex]
                            * elementUniversal.globalData.Wc_2D[integrationPointIndex / elementUniversal.globalData.nPc]
                            * elementUniversal.globalData.Wc_2D[integrationPointIndex % elementUniversal.globalData.nPc]
                            * J_i_matrices[integrationPointIndex].Determinant();

                        C_l_matrices[gridElementIndex] += C_i_matrices[integrationPointIndex]
                            * elementUniversal.globalData.Wc_2D[integrationPointIndex / elementUniversal.globalData.nPc]
                            * elementUniversal.globalData.Wc_2D[integrationPointIndex % elementUniversal.globalData.nPc]
                            * J_i_matrices[integrationPointIndex].Determinant();
                    }
                }


                // uwzględniam warunki brzegowe w macierzy H_l i wektorze P_l dla każdego elementu
                for (int gridElementIndex = 0; gridElementIndex < grid.nE; gridElementIndex++)
                {
                    Element element = grid.Elements[gridElementIndex];
                    var sidesWithBcIndexes = new List<int>();

                    // sprawdzam, czy są boki elementu z bc = 1, jeżeli dany bok ma bc to dodaję jego indeks do listy indeksów
                    for (int sideIndex = 0; sideIndex < element.SidesWithBc.Count; sideIndex++)
                    {
                        if (element.SidesWithBc[sideIndex] == true)
                        {
                            sidesWithBcIndexes.Add(sideIndex);
                        }
                    }

                    if (sidesWithBcIndexes.Count > 0)
                    {
                        // macierz H_bc_l dla danego elementu, którą będziemy dodawać do macierzy H_l tego elementu
                        var H_bc_l_matrix = Matrix<double>.Build.Dense(numberOfNodesInElement, numberOfNodesInElement);

                        // wektor P_l dla danego elementu
                        var P_l_vector = Matrix<double>.Build.Dense(numberOfNodesInElement, 1);

                        for (int sideWithBcIndex = 0; sideWithBcIndex < sidesWithBcIndexes.Count; sideWithBcIndex++)
                        {
                            int nPc = elementUniversal.globalData.nPc;

                            // lista macierzy H_bc_i dla każdego punktu całkowania leżącego na danym boku
                            var H_bc_i_matrices = new List<Matrix<double>>(nPc);

                            // lista wektorów P_i dla każdego punktu całkowania leżącego na danym boku
                            var P_i_vectors = new List<Matrix<double>>(nPc);


                            // dla każdego punktu całkowania tworzę macierz H_bc_i i wektor P_i
                            for (int integrationPointIndex = 0; integrationPointIndex < nPc; integrationPointIndex++)
                            {
                                double ksi = elementUniversal.globalData.Pc_1D_Ksi[integrationPointIndex + nPc * sidesWithBcIndexes[sideWithBcIndex]];
                                double eta = elementUniversal.globalData.Pc_1D_Eta[integrationPointIndex + nPc * sidesWithBcIndexes[sideWithBcIndex]];

                                var N_i = Matrix<double>.Build.Dense(numberOfNodesInElement, 1);

                                N_i[0, 0] = 0.25 * (1 - ksi) * (1 - eta);
                                N_i[1, 0] = 0.25 * (1 + ksi) * (1 - eta);
                                N_i[2, 0] = 0.25 * (1 + ksi) * (1 + eta);
                                N_i[3, 0] = 0.25 * (1 - ksi) * (1 + eta);

                                Matrix<double> H_bc_i_matrix = alfa * (N_i * N_i.Transpose());
                                Matrix<double> P_i_vector = alfa * N_i * t_alfa;

                                H_bc_i_matrices.Add(H_bc_i_matrix);
                                P_i_vectors.Add(P_i_vector);

                                //Console.WriteLine(string.Format("element nr: {0}, bok nr: {1}, wektor P_i dla pc nr: {2}", gridElementIndex, sidesWithBcIndexes[sideWithBcIndex], integrationPointIndex));
                                //Console.WriteLine(P_i_vector.ToString());
                            }

                            // sumuję macierze H_bc_i i wektory P_i dla każdego z punktów całkowania i dodaję je do odpowiednich list
                            for (int integrationPointIndex = 0; integrationPointIndex < H_bc_i_matrices.Capacity; integrationPointIndex++)
                            {
                                double x_1 = element.Nodes[sidesWithBcIndexes[sideWithBcIndex]].x;
                                double x_2 = element.Nodes[(sidesWithBcIndexes[sideWithBcIndex] + 1) % numberOfNodesInElement].x;
                                double y_1 = element.Nodes[sidesWithBcIndexes[sideWithBcIndex]].y;
                                double y_2 = element.Nodes[(sidesWithBcIndexes[sideWithBcIndex] + 1) % numberOfNodesInElement].y;

                                double deltaX = Math.Sqrt(Math.Pow(x_1 - x_2, 2) + Math.Pow(y_1 - y_2, 2));

                                H_bc_l_matrix += H_bc_i_matrices[integrationPointIndex] * elementUniversal.globalData.Wc_2D[integrationPointIndex] *
                                    (deltaX / 2);

                                P_l_vector += P_i_vectors[integrationPointIndex] * elementUniversal.globalData.Wc_2D[integrationPointIndex] *
                                    (deltaX / 2);
                            }
                        }

                        H_bc_l_matrices[gridElementIndex] += H_bc_l_matrix;
                        P_l_vectors[gridElementIndex] += P_l_vector;
                    }
                }

                //Console.WriteLine("\n\n///////////////////////////////////////////////////////////////");
                //Console.WriteLine("// Macierze H_l, C_l, wektory P_l dla odpowiednich elementów //");
                //Console.WriteLine("///////////////////////////////////////////////////////////////\n\n");

                for (int gridElementIndex = 0; gridElementIndex < grid.nE; gridElementIndex++)
                {
                    //Console.WriteLine("/////////////\nElement nr: " + (gridElementIndex + 1));

                    //PrintHGCP(H_l_matrices[gridElementIndex], (H_l_matrices[gridElementIndex] + H_bc_l_matrices[gridElementIndex]),
                    //    C_l_matrices[gridElementIndex], P_l_vectors[gridElementIndex]);
                }


                ///////////////////////////////////////////////////
                ///// Agregacja:                              /////
                /////  - macierzy H_l z war.brzeg. do H_g,    /////
                /////  - macierzy C_l do C_g,                 /////
                /////  - wektorów P_l do P_g.                 /////
                ///////////////////////////////////////////////////

                // macierz H globalna
                var H_g_matrix = Matrix<double>.Build.Dense(grid.nN, grid.nN);

                // macierz C globalna
                var C_g_matrix = Matrix<double>.Build.Dense(grid.nN, grid.nN);

                // wektor P globalny
                var P_g_vector = Matrix<double>.Build.Dense(grid.nN, 1);


                // agregacja macierzy H_g, C_g i wektora P_g
                for (int gridElementIndex = 0; gridElementIndex < grid.nE; gridElementIndex++)
                {
                    Element element = grid.Elements[gridElementIndex];

                    var IDs = new List<int>();

                    for (int nodeIndex = 0; nodeIndex < numberOfNodesInElement; nodeIndex++)
                    {
                        //Console.WriteLine(string.Format("element nr: {0}, nodeIndexL: {1}, nodeIndexG: {2}", gridElementIndex, nodeIndex, element.Nodes[nodeIndex].ID));
                        IDs.Add(element.Nodes[nodeIndex].ID - 1);
                    }

                    for (int i = 0; i < IDs.Count; i++)
                    {
                        P_g_vector[IDs[i], 0] += P_l_vectors[gridElementIndex][i, 0];

                        for (int j = 0; j < IDs.Count; j++)
                        {
                            H_g_matrix[IDs[i], IDs[j]] += H_l_matrices[gridElementIndex][i, j] + H_bc_l_matrices[gridElementIndex][i, j];
                            C_g_matrix[IDs[i], IDs[j]] += C_l_matrices[gridElementIndex][i, j];
                        }
                    }
                }

                Console.WriteLine("\n\n///////////////////////////////////////////////////////////////");
                Console.WriteLine("//  Zagregowane: Macierze H_g, macierz C_g, wektor P_g       //");
                Console.WriteLine("///////////////////////////////////////////////////////////////\n\n");

                Console.WriteLine("Macierz H globalna");
                Console.WriteLine(H_g_matrix.ToString());

                Console.WriteLine("Macierz C globalna");
                Console.WriteLine(C_g_matrix.ToString());

                Console.WriteLine("Wektor P globalny");
                Console.WriteLine(P_g_vector.ToString());


                /////////////////////////////////////////////////////////
                ///// Rozwiązuję układ równań:                      /////
                /////  - otrzymuję {t_1},                           /////
                /////  - podstawiam {t_0} wartościami {t_1},        /////
                /////  - powtarzam algorytm tyle, ile jest kroków.  /////
                /////////////////////////////////////////////////////////            

                // tworzę wektory {t_0} i {t_1}
                var T_0_vector = Matrix<double>.Build.Dense(grid.nN, 1);
                var T_1_vector = Matrix<double>.Build.Dense(grid.nN, 1);

                Console.WriteLine("\n\n///////////////////////////////////////////////////////////////");
                Console.WriteLine("//////////////////// Rozwiązanie: /////////////////////////////");
                Console.WriteLine("///////////////////////////////////////////////////////////////\n\n");

                for (int iteration = 0; iteration < nIterations; iteration++)
                {
                    // węzły siatki nie dodają się po kolei względem wartości ID
                    for (int nodeIndex = 0; nodeIndex < grid.Nodes.Count; nodeIndex++)
                    {
                        T_0_vector[grid.Nodes[nodeIndex].ID - 1, 0] = grid.Nodes[nodeIndex].t;
                    }

                    var H_g_1_matrix = H_g_matrix + C_g_matrix.Divide(dTau);
                    var P_g_1_vector = P_g_vector + C_g_matrix.Divide(dTau) * T_0_vector;

                    T_1_vector = H_g_1_matrix.Inverse() * P_g_1_vector;

                    // szukam temperatury minimalnej i maksymalnej
                    double minTemp = T_1_vector[0, 0];
                    double maxTemp = T_1_vector[0, 0];

                    for (int i = 0; i < T_1_vector.RowCount; i++)
                    {
                        if (T_1_vector[i, 0] >= maxTemp) maxTemp = T_1_vector[i, 0];
                        if (T_1_vector[i, 0] < minTemp) minTemp = T_1_vector[i, 0];
                    }


                    Console.WriteLine("\n\n///////////////////////////////////////////////");
                    Console.WriteLine("iteracja nr: " + iteration + "\n");

                    Console.WriteLine("\nMacierz [H] = [H] + [C]/dTau");
                    Console.WriteLine(H_g_1_matrix.ToString());

                    Console.WriteLine("\nWektor {p} = {p} + {[C]/dTau} * {t_0}");
                    Console.WriteLine(P_g_1_vector.ToString());

                    Console.WriteLine(string.Format("\nCzas[s] = {0}, MinTemp[C] = {1}, MaxTemp[C] = {2}\n\n", (iteration + 1) * dTau, minTemp, maxTemp));

                    for (int nodeIndex = 0; nodeIndex < grid.Nodes.Count; nodeIndex++)
                    {
                        grid.Nodes[nodeIndex].t = T_1_vector[grid.Nodes[nodeIndex].ID - 1, 0];
                    }

                    if (iteration == numberOfIterationsToPrint) break;
                }



                Console.ReadKey();
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                Console.ReadKey();
            }
        }

        static void PrintHGCP(Matrix<double> H, Matrix<double> H_Hbc, Matrix<double> C, Matrix<double> P)
        {
            Console.WriteLine("\tmacierz H_l:\n");
            Console.WriteLine(H.ToString());
            Console.WriteLine("\tmacierz H_l + H_bc_l:\n");
            Console.WriteLine(H_Hbc.ToString());
            Console.WriteLine("\tmacierz C_l:\n");
            Console.WriteLine(C.ToString());
            Console.WriteLine("\twektor P_l:\n");
            Console.WriteLine(P.ToString());
            Console.WriteLine("/////////////////////////////////////////////////\n");
        }
    }
}
