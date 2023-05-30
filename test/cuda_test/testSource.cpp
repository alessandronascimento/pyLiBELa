#include "cudaTest.cuh"

void cudaWrapper(){
  execution();
}

// int main() {
//   execution();
//   return 0;
// }

#include "boost/python.hpp"

BOOST_PYTHON_MODULE(cudaTest) {
  using namespace boost::python;
  def("cudaTest", cudaWrapper);
}

