/* Matrix addition example

  Copyright (c) 2024 Debajyoti Debnath                                         
                                                                           
  Licensed under the Apache License, Version 2.0 (the "License");          
  you may not use this file except in compliance with the License.         
  You may obtain a copy of the License at                                  
                                                                           
      http://www.apache.org/licenses/LICENSE-2.0                           
                                                                           
  Unless required by applicable law or agreed to in writing, software      
  distributed under the License is distributed on an "AS IS" BASIS,        
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
  See the License for the specific language governing permissions and      
  limitations under the License.                                           
*/

#include <stdio.h>
#include "matrix.h"

int main() {
    size_t seed = 42;
    int nrow = 10, ncol = 5;
    Matrix a = mat_create(nrow, ncol);
    mat_fill_random(&a, seed);
    Matrix b = mat_copy(&a);
    mat_scale(&b, 2.34);
    
    printf("\n************Matrix addition***************\n");
    printf("A: \n");
    mat_print(&a);
    printf(" +\n B: \n");
    mat_print(&b);
    printf(" = \n");
    
    // Addition
    mat_add(&a, &b);
    mat_print(&a);
    printf("\n******************************************\n");

    mat_destroy(&a);
    mat_destroy(&b);
}
