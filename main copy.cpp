
#include <cstdio>
#include <random>
#include <string>
#include <iostream>
#include <vector>
#include <numeric>
#include <float.h>
#include <math.h>

#include <GridCut/GridGraph_2D_4C.h>

void basic()
{

  // api定义见库目录下readme

  typedef GridGraph_2D_4C<int, int, int> Grid;

  Grid *grid = new Grid(2, 2);

  grid->set_terminal_cap(grid->node_id(0, 0), 10, 0);
  grid->set_terminal_cap(grid->node_id(0, 1), 0, 10);
  grid->set_terminal_cap(grid->node_id(1, 0), 10, 0);
  grid->set_terminal_cap(grid->node_id(1, 1), 0, 10);

  grid->set_neighbor_cap(grid->node_id(0, 0), +1, 0, 1);
  grid->set_neighbor_cap(grid->node_id(0, 0), 0, +1, 2);

  grid->set_neighbor_cap(grid->node_id(1, 0), -1, 0, 5);
  grid->set_neighbor_cap(grid->node_id(1, 0), 0, +1, 3);

  grid->set_neighbor_cap(grid->node_id(0, 1), 0, -1, 2);
  grid->set_neighbor_cap(grid->node_id(0, 1), +1, 0, 1);

  grid->set_neighbor_cap(grid->node_id(1, 1), 0, -1, 3);
  grid->set_neighbor_cap(grid->node_id(1, 1), -1, 0, 1);

  grid->compute_maxflow();

  printf("Min-cut partition:\n");

  for (int x = 0; x < 2; x++)
  {
    for (int y = 0; y < 2; y++)
    {
      printf("%c", (grid->get_segment(grid->node_id(x, y)) == 0) ? 'S' : 'T');
    }
    printf("\n");
  }

  delete grid;
}

void testPerformance(int size)
{
  std::default_random_engine e;
  std::uniform_int_distribution<int> u(0, 100);
  e.seed(time(0));

  typedef GridGraph_2D_4C<int, int, int> Grid;

  Grid *grid = new Grid(size, size);

  for (int x = 0; x < size - 1; x++)
  {
    for (int y = 0; y < size - 1; y++)
    {
      grid->set_terminal_cap(grid->node_id(x, y), u(e), u(e));
      grid->set_neighbor_cap(grid->node_id(x, y), +1, 0, u(e));
      grid->set_neighbor_cap(grid->node_id(x, y), 0, +1, u(e));
    }
  }

  grid->compute_maxflow();
  printf("Min-cut partition\n");
  for (int x = 0; x < size; x++)
  {
    for (int y = 0; y < size; y++)
    {
      printf("%c", (grid->get_segment(grid->node_id(x, y)) == 0) ? 'S' : 'T');
    }
    printf("\n");
  }
  delete grid;
}

class Polygon
{
public:
  std::string name;
  std::string label;
  double height[3];
  Polygon() {}

  Polygon(std::string name, std::string label, double h1, double h2, double h3)
  {
    this->name = name;
    this->label = label;
    this->height[0] = h1;
    this->height[1] = h2;
    this->height[2] = h3;
  }

  ~Polygon()
  {
  }
};

struct Label
{
  std::string label;
  double planeHeight;
  std::vector<double> pointCollection;
};

// init fake roof
Polygon pList[9]{
    Polygon("poly00", "p1", 10.2, 10.12, 9.9),
    Polygon("poly01", "p2", 15.22, 15.16, 14.87),
    Polygon("poly02", "p3", 10.1, 10.01, 9.95),
    Polygon("poly10", "p4", 10.02, 10.15, 9.87),
    Polygon("poly11", "p4", 5.2, 5.1, 4.9),
    Polygon("poly12", "p3", 10.02, 10.15, 9.87),
    Polygon("poly20", "p2", 10.2, 10.1, 9.9),
    Polygon("poly21", "p4", 10.1, 10.01, 9.95),
    Polygon("poly22", "p2", 10.2, 10.1, 9.9),
};

Label f[4];
void countPointInLabel()
{
  f[0].label == "p1";
  f[0].label == "p2";
  f[0].label == "p3";
  f[0].label == "p4";
  for (int i = 0; i < 9; i++)
  {
    if (pList[i].label == "p1")
    {
      f[0].pointCollection.push_back(pList[i].height[0]);
      f[0].pointCollection.push_back(pList[i].height[1]);
      f[0].pointCollection.push_back(pList[i].height[2]);
    }
    else if (pList[i].label == "p2")
    {
      f[1].pointCollection.push_back(pList[i].height[0]);
      f[1].pointCollection.push_back(pList[i].height[1]);
      f[1].pointCollection.push_back(pList[i].height[2]);
    }
    else if (pList[i].label == "p3")
    {
      f[2].pointCollection.push_back(pList[i].height[0]);
      f[2].pointCollection.push_back(pList[i].height[1]);
      f[2].pointCollection.push_back(pList[i].height[2]);
    }
    else if (pList[i].label == "p4")
    {
      f[3].pointCollection.push_back(pList[i].height[0]);
      f[3].pointCollection.push_back(pList[i].height[1]);
      f[3].pointCollection.push_back(pList[i].height[2]);
    }
  }
}

void energyMinimization()
{

  /*
 1. Start with an arbitrary labeling f

 2. For each label a in L
 2.1. Find f' = argminE(f0) among f0 within one expansion of f
 2.2  If E(f') < E(f), set f = f' and goto 3

 3. Return f  //no more f' founded
 */

  // calc E(f)

  countPointInLabel();

  for (int l = 0; l < 4; l++)
  {
    // label=>p(l)
    // for label p
    // get all points marked in p
    std::vector<double> v = f[l].pointCollection;
    double planeHeight = std::accumulate(std::begin(v), std::end(v), 0.0) / v.size(); // meanheight
    f[l].planeHeight = planeHeight;

    // calc Dp
    double Dp = 0;
    for (int pi = 0; pi < v.size(); pi++)
    {
      Dp += (v[pi] - planeHeight) * (v[pi] - planeHeight);
    }
    // calc Vpq
    /*
    00-01-02
    10-11-12
    20-21-22
    */
    double Vpq = 0;
    for (int row = 0; row < 2; row++)
    {
      for (int col = 0; col < 2; col++)
      {
        if (pList[row * 3 + col].label != pList[row * 3 + col + 1].label)
          Vpq += 1;
        if (pList[row * 3 + col].label != pList[(row + 1) * 3 + col].label)
          Vpq += 1;
      }
    }
    if (
        pList[2].label != pList[5].label ||
        pList[5].label != pList[8].label ||
        pList[6].label != pList[7].label ||
        pList[7].label != pList[8].label)
    {
      Vpq += 1;
    }

    std::cout << "label:" << f[l].label
              << " planeHeight:" << planeHeight
              << " D(p):" << Dp
              << " V(pq):" << Vpq
              << std::endl;

    // construct alpha-beta swap graph

    typedef GridGraph_2D_4C<int, int, int> Grid;

    Grid *grid = new Grid(3, 3);

    Polygon p = pList[0];
    Polygon q = pList[1];

    double tpl, tql, tpnl, tqnl; // Dp when l is assigned to p

    tpl = std::pow(f[l].planeHeight - p.height[0], 2) + std::pow(f[l].planeHeight - p.height[1], 2) + std::pow(f[l].planeHeight - p.height[2], 2);
    tql = std::pow(f[l].planeHeight - q.height[0], 2) + std::pow(f[l].planeHeight - q.height[1], 2) + std::pow(f[l].planeHeight - q.height[2], 2);

    if (f[l].label == p.label)
    {
      double tpnl = DBL_MAX;
    }
    else
    {
      for (int label = 0; label < 4; label++)
      {
        if (f[label].label == p.label)
        {
          tpnl = std::pow(f[label].planeHeight - p.height[0], 2) + std::pow(f[label].planeHeight - p.height[1], 2) + std::pow(f[label].planeHeight - p.height[2], 2);
        }
      }
    }

    if (f[l].label == q.label)
    {
      double tqnl = DBL_MAX;
    }
    else
    {
      for (int label = 0; label < 4; label++)
      {
        if (f[label].label == q.label)
        {
          tpnl = std::pow(f[label].planeHeight - q.height[0], 2) + std::pow(f[label].planeHeight - q.height[1], 2) + std::pow(f[label].planeHeight - q.height[2], 2);
        }
      }
    }

    double epa,eaq,tanl;
    if(p.label!=q.label){
      epa=p.label==f[l].label?0:1;
      eaq=q.label==f[l].label?0:1;
      tanl=1;
    }else{
      epa=DBL_MAX;
      eaq=DBL_MAX;
      tanl=1;
    }
   

  }
  
}


void simplRoofCut()
{
}

int main(int argc, char **argv)
{

  energyMinimization();
  return 0;
}