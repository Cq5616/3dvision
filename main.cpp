

#include <string>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <random>
#include <memory>

#include <AlphaExpansion/AlphaExpansion_2D_4C.h>

class Polygon
{
public:
  std::string name;
  double groundTruth;
  std::vector<double> sampleHeight;
  Polygon() {}
  Polygon(std::string name, double groundTruth) : name(name), groundTruth(groundTruth) {}

  void generateSamples(int num, double var)
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    for (int i = 0; i < num; i++)
    {
      sampleHeight.push_back(groundTruth + U(gen) * var);
    }
  }

  ~Polygon()
  {
  }
};

const int GraphWidth = 3, GraphHeight = 3, numLabels = 2;
int numSamples = 10;
struct LabelList
{
  std::vector<double> sampleHeight[numLabels];
  double planeHeight[numLabels];
};

void estimatePlaneHeight(std::vector<Polygon> polyList, int *labeling, std::shared_ptr<LabelList> lList)
{

  for (int i = 0; i < polyList.size(); i++)
  {
    int label = labeling[i];
    for (int s = 0; s < polyList[i].sampleHeight.size(); s++)
    {
      lList->sampleHeight[label].push_back(polyList[i].sampleHeight[s]);
    }
  }

  for (int i = 0; i < numLabels; i++)
  {
    if (lList->sampleHeight[i].size() != 0)
    {
      lList->planeHeight[i] = std::accumulate(lList->sampleHeight[i].begin(), lList->sampleHeight[i].end(), 0.0) / lList->sampleHeight[i].size();
    }
    else{
      lList->planeHeight[i]=0.0;
    }
    std::cout << "label" << i << ": sample num:" << lList->sampleHeight[i].size() << " height:" << lList->planeHeight[i] << std::endl;
  }
}

double smoothFn(int pix1, int pix2, int label1, int label2)
{
  if (label1 == label2)
    return 0.0;
  return 1.5;
}

double dataCostFn(Polygon poly,int label,std::shared_ptr<LabelList> lList){
  double datacost=0.0;
  for(int i=0;i<poly.sampleHeight.size();i++){
    datacost+=std::abs(poly.sampleHeight[i]-lList->planeHeight[label]);
  }
  datacost/=poly.sampleHeight.size();
  return datacost;
}

void energyMinimization(std::vector<Polygon> polyList, int *labeling)
{
  auto lList= std::make_shared<LabelList>();
  
  std::cout << "----estimatePlaneHeight----" << std::endl;
  // init plane by current labeling
  estimatePlaneHeight(polyList,labeling, lList);

  // int label [0,1,2,3] => plane0,plane1,plane2,plane3
  typedef AlphaExpansion_2D_4C<int, double, double> Expansion;

  // dataCosts
  /*
  We expect that data costs are provided by the array of n*k elements, i.e. the array contains k values for each pixel.
    The values are ordered in the way, that for the pixel index p_i and the label l, the data cost value is
    on the array index p_i*k+l.

      l_0(p_0) | l_1(p_0) | ... | l_{k-1}(p_0) | l_0(p_1) | l_1(p_1) | ... | l_{k-1}(p_1) | ... | ... | l_{k-1}(p_{n-1})

  */
  
  double dataCosts[GraphWidth * GraphHeight * numLabels];
  for (int row = 0; row < GraphWidth; row++)
  {
    for (int col = 0; col < GraphHeight; col++)
    {
      for (int l = 0; l < numLabels; l++)
      {
        int polyIdx = (row * GraphWidth + col) ;
      
        dataCosts[polyIdx* numLabels + l] = dataCostFn(polyList[polyIdx],l,lList);

        //std::cout << "when set l" << l << " (planeheight:" << lList->planeHeight[l] << ") to polygon" << row << col << " cost:" << dataCosts[polyIdx + l] << std::endl;
      }
    }
  }

  // Graph Cuts Energy Minimization
  Expansion *expansion = new Expansion(GraphWidth, GraphHeight, numLabels, dataCosts, &smoothFn);
  expansion->set_labeling(labeling);

  expansion->perform_random();
  std::cout << "----after Energy Minimization----" << std::endl;
  std::cout << "energy:" << expansion->get_energy() << std::endl;
  labeling = expansion->get_labeling();

}

void iterator(std::vector<Polygon> polyList)
{
  int steps = 10;
  int *labeling = new int[GraphWidth * GraphHeight];
  //init random labeling
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> U(0,1.0);
  for(int i=0;i<GraphWidth*GraphHeight;i++){
    labeling[i]=std::floor( U(gen)/(1.0/numLabels) );
  }

  //it
  for (int i = 0; i < steps; i++)
  {
    std::cout << "===iterator step:" << i << "===" << std::endl;
    energyMinimization(polyList, labeling);
    std::cout << "labeling: ";
    for (int j = 0; j < GraphWidth * GraphHeight; j++)
    {
      std::cout << labeling[j] << " ";
    }
    std::cout << std::endl;
  }

}

int main(int argc, char **argv)
{
  std::vector<Polygon> polyList;
  // init fake data
  std::cout<<"ground truth: "<<std::endl;
  for (int row = 0; row < GraphHeight; row++)
  {
    for (int col = 0; col < GraphWidth; col++)
    {
      Polygon p("poly" + row + col, row > col ? 100 : 50); // ground truth 对角阵列
      p.generateSamples(numSamples, 2.0);
      polyList.push_back(p);

      std::cout<<p.groundTruth<<" ";
    }
  }
  std::cout<<std::endl;
  iterator(polyList);
  
  return 0;
}