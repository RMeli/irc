/* 
    Author: hauptmech <hauptmech@gmail.com>, Nov 2013
    Modified by Rocco Meli, Feb 2018

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

EIGEN_STRONG_INLINE Matrix(std::initializer_list<_Scalar> initlist) : Base()
{
  
  // Check template parameters
  Base::_check_template_params();
  
  // Get std::initializer_list size
  std::size_t size{initlist.size()};
  
  // Resize array, if empty
  if (base().size() == 0) {
    // TODO: Assumes column vector
    this->resize(size, 1);
  }
  else{
    // Check size compatibility between matrix and INITLIST
    eigen_assert(size == base().size());
  }
  
  // Fill array from INITLIST
  std::size_t i{0};
  for(const auto& x : initlist){
    this->operator[](i) = x;
    i++;
  }
}

EIGEN_STRONG_INLINE Matrix(std::initializer_list<std::initializer_list<_Scalar>> initlist) : Base()
{
  // Check template parameters
  Base::_check_template_params();
  
  std::size_t rows{initlist.size()};
  std::size_t cols{initlist.begin()->size()};
  
  // Resize array, if empty
  if (base().size() == 0) {
    this->resize(rows, cols);
  }
  else{
    // Check size compatibility between matrix and INITLIST
    eigen_assert(rows == _Rows);
    eigen_assert(cols == _Cols);
  }
  
  // Fill array from INITLIST
  std::size_t i{0};
  for(const auto& x : initlist){
    std::size_t j{0};
    for(const auto& y : x){
      this->operator()(i,j) = y;
      j++;
    }
    i++;
  }
}
