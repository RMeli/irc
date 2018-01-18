/* 
    Author: hauptmech <hauptmech@gmail.com>, Nov 2013 

    This is free and unencumbered software released into the public domain.

    Anyone is free to copy, modify, publish, use, compile, sell, or
    distribute this software, either in source code form or as a compiled
    binary, for any purpose, commercial or non-commercial, and by any
    means.

    In jurisdictions that recognize copyright laws, the author or authors
    of this software dedicate any and all copyright interest in the
    software to the public domain. We make this dedication for the benefit
    of the public at large and to the detriment of our heirs and
    successors. We intend this dedication to be an overt act of
    relinquishment in perpetuity of all present and future rights to this
    software under copyright law.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
    OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

    For more information, please refer to <http://unlicense.org/>
*/

#include <initializer_list>


/** Constructs an array from a C++11 initializer_list

  ArrayXXd a {1,2,3,4} //4x1 array
  ArrayXd b {6,5,8,6} // 4 element array
  Array3d c {1,2,3}

  The constructor will have an assertion failure if a fixed size
  array given a initializer list of the wrong size.

 **/
EIGEN_STRONG_INLINE Array(std::initializer_list<double > initlist) : Base()
{

  Base::_check_template_params();
  int size = initlist.size();

  if (base().size() == 0) { //Empty array, gotta resize it
    this->resize(size,1 );
  }

  int initializer_list_size = initlist.size();
  int Array_size = base().size();
  eigen_assert(initializer_list_size == Array_size);

  if (static_cast<int>(initlist.size()) == base().size()){
    int i = 0;
    for (auto x: initlist){
       coeffRef(i)=Scalar(x);
       i++;
    }
  }
}

/** Constructs an array from a C++11 initializer_list

  ArrayXXd a {1,2,3,4} //4x1 array
  ArrayXd b {6,5,8,6} // 4 element array
  Array3d c {1,2,3}

  The constructor will have an assertion failure if a fixed size
  array given a initializer list of the wrong size.

 **/
EIGEN_STRONG_INLINE Array(std::initializer_list<std::initializer_list<double> > initlist) : Base()
{
  Base::_check_template_params();
  int rows = initlist.size();
  int cols = initlist.begin()->size();

  if (base().size() == 0) { //Empty array, gotta resize it
    this->resize(rows, cols);
  }

  int initializer_list_size = rows*cols;
  int Array_size = base().size();
  eigen_assert(initializer_list_size == Array_size);

  if (rows*cols == base().size()){
    int i = 0;
    for (auto x: initlist){
      int j=0;
      for (auto y: x){
        coeffRef(rows*j+i)=Scalar(y);
        j++;
      }
      i++;
    }
  }
}
