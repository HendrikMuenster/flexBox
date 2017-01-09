#include <iostream>
#include <omp.h>

#include "flexBox.h"

using namespace std;

int main()
{
  cout << "Hello World !!"  << endl;


  //check flexbox
  flexBox<float,std::vector<float>> mainObject;

  //check openmp link if existing
  int numThreads, id;
  #pragma omp parallel private(numThreads, id)
  {
    id = omp_get_thread_num();

    #pragma omp critical
    cout << "Thread: " << id << endl;

    if(id == 0)
    {
      numThreads = omp_get_num_threads();
      #pragma omp critical
      cout << "Num threads: " << numThreads << endl;
    }
  }

  return 0;
}
