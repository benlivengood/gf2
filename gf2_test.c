/* Copyright (c) 2013, Ben Livengood */
/* All rights reserved. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions are met: */
/* * Redistributions of source code must retain the above copyright */
/*       notice, this list of conditions and the following disclaimer. */
/*       * Redistributions in binary form must reproduce the above copyright */
/*       notice, this list of conditions and the following disclaimer in the */
/*       documentation and/or other materials provided with the distribution. */
/*       * Neither the name of the copyright owner nor the */
/*       names of its contributors may be used to endorse or promote products */
/*       derived from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND */
/* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED */
/* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE */
/* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY */
/* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES */
/* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT */
/* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include "gf2.h"

#define TEST_ITER 1000

gf2 test_field;

gf2_matrix identity2x2={rows:2,columns:2,field:&test_field,
			element:{0x01,0x00,
				 0x00,0x01}}; /* 2x2 identity matrix */

gf2_matrix test4x2={rows:4,columns:2,field:&test_field,
		    element:{0x1,0x2,0x3,0x4,
			     0x5,0x6,0x7,0x8}
}; /* 4x2 matrix laid out in transposed form */

gf2_matrix test2x4={rows:2,columns:2,field:&test_field,
		    element:{0x9,0xa,
			     0xb,0xc,
			     0xd,0xe,
			     0xf,0x10}
}; /* 2x4 matrix laid out in transposed form */


gf2_matrix testcolumn={rows:4,columns:1,field:&test_field,
		       element:{0x20,0x21,0x22,0x23}}; /* 4x1 matrix */

gf2_matrix testrow={rows:1,columns:4,field:&test_field,
		    element:{0x30,0x31,0x32,0x33}}; /* 1x4 matrix */



int test_gf2_matrix()
{
  int result;
  printf("\nTesting gf2_matrix functions\n");


  if ((result=gf2_init(&test_field,0x1d)) !=0 ) {
    printf("test_gf2_matrix: gf2_init(&test)field,0x1d) failed on element %d\n",result+1000);
    return 1;
  }

  gf2_matrix *r=0;

  r=gf2_matrix_mul(&test4x2,&identity2x2,r);
  if (!r) {
    printf("test_gf2_matrix: gf2_matrix_mul returned zero\n");
    return 2;
  }
  printf("r=\n");
  gf2_matrix_print(r);



  return 0;
}


int gcd(int a,int b)
{
  int c;
  while(b) {
    c=a%b;
    a=b;
    b=c;
  }
  return a;
}

int totient(int a)
{
  int i,t=0;
  for(i=1;i<a;i++) {
    if (gcd(i,a)==1) t++;
  }
  return t;
}

int primitive_polynomials(int p,int m)
{
  int i,d=1;
  for(i=0;i<m;i++) d*=p;
  d--;
  return totient(d)/m;
}


int main(int argc,char **argv)
{
  gf2 field[1]; /* Field to test */
  int i,result,r,ir,j,k;
  printf("Testing gf2...\n");
  
  for(i=1;i<12;i++) {
    printf("Primitive polynomials in GF(2^%d): %d\n",i,
	   primitive_polynomials(2,i)); 
  }

  /* Exhaustive test of all 8th degree polynomials */  
  printf("Testing all 8th degree polynomials in GF(2)\n");
  r=0;
  ir=0;
  for(i=0;i<256;i++) {
    result=gf2_init(field,i); /* Test polynomial i */
    if (result==0) {
      ir++;
      printf("1%.2x is irreducible.\n",i);
      /* Test symmetry of multiplication and equality of both methods */
      for(j=0;j<256;j++) {
	for(k=0;k<256;k++) {
	  if (gf2_mul1(field,j,k)!=gf2_mul1(field,k,j) ||
	      gf2_mul2(field,j,k)!=gf2_mul1(field,k,j)) {
	    printf("FAIL: %.2x * %.2x does not equal %.2x * %.2x\n",i,j,j,i);
	  }
	}
      }
      /* Test logarithms and exponents */
      for(j=1;j<255;j++) {
	if (gf2_exp(field,2,gf2_log(field,j))!=
	    gf2_log(field,gf2_exp(field,2,j))) {
	  printf("FAIL: exp(log(%.2x)) (%.2x) != log(exp(%.2x)) (%.2x)\n",j,j,
		 gf2_exp(field,2,gf2_log(field,j)),
		 gf2_log(field,gf2_exp(field,2,j)));
	}
      }
      /* Test edge case logarithms and exponents */

      for(j=0;j<256;j++) {
	if (gf2_mul1(field,0,j)!=0 ||
	    gf2_mul2(field,j,0)!=0) {
	  printf("FAIL: Multiplication by zero doesn't yield zero %d\n",j);
	}
      }

      if (gf2_exp(field,2,gf2_log(field,0xff))!=0xff)
	  printf("FAIL: exp(log(ff)) != ff\n");
      if (gf2_log(field,gf2_exp(field,2,0xff))!=0x00)
	  printf("FAIL: loglog(ff)) != ff\n");

      /* Multiply all nonzero elements together which will produce 1. This is
      ** because every nonzero element has an inverse, and except
      ** for 1 each inverse is not itself. So every element will be multiplied
      ** by its inverse only once. */
      j=1;
      for(k=1;k<256;k++) j=gf2_mul2(field,j,k);
      if (j!=1) printf("FAIL: Product of all nonzero elements != 01\n");

      /* Make sure all inverses times themselves equal 1 */
      for(k=1;k<256;k++) {
	if (gf2_mul2(field,k,gf2_inv(field,k))!=1) {
	  printf("FAIL: %.2x * %.2x != 01\n",k,gf2_inv(field,k));
	}
      }

      
      /*      Make sure exponents are working by checking Fermat's theorem */
      for(k=1;k<256;k++) {
	for(j=0;j<256;j++) {
	  if (gf2_mul2(field,gf2_exp(field,k,j),gf2_exp(field,k,256-j))!=k)
	    printf("FAIL: %.2x^%d * %.2x^%d != %.2x\n",k,j,k,256-j,k);
	}
      }
    }
    else {
      r++;
    }
  }
  if (r+ir != 256) printf("FAIL: Only %d polynomials tested\n",r+ir);
  printf("Reducible polynomials:   %.3d\n",r);
  printf("Irreducible polynomials: %.3d\n",ir);

  if (primitive_polynomials(2,8)!=ir) {
    printf("FAIL: Should be %d primitive polynomials over GF(2^8)\n",
	   primitive_polynomials(2,8)); 
  }

  for(i=0;i<256;i++) if (gf2_init(field,i)==0) break;
  if (i==256) {
    printf("FAIL: Couldn't find an irreducible polynomial\n");
    return -1;
  }


  if (test_gf2_matrix() != 0) {
    printf("FAIL: test_gf2_matrix returned nonzero\n");
    return 100;
  }

  
  gf2_matrix *vand,*generator;
  
  int codesize=255,datasize=128;

  srand(time(0));

  printf("\n");

  vand=gf2_matrix_alloc(codesize,datasize,field);
  if (vand==0) {
    printf("FAIL: Out of memory: vand\n");
    return 1;
  }

  for(i=0;i<datasize;i++) {
    for(j=0;j<codesize;j++) {
      vand->element[j+i*codesize]=gf2_exp(field,j,i);
    }
  }

  generator=gf2_gaussian_elimination(vand,0);
  if (generator==0) {
    printf("FAIL: gaussian elimination failed\n");
    return 1;
  }

  




  gf2_matmul *mm;
  gf2_mulmat *mm2;
  
  gf2_matrix *a,*b,*c,*d,*e,*f,*g,*h;

  srand(time(0));

  a = gf2_matrix_alloc(16,16,field);
  for(i=0;i<16*16;i++) a->element[i]=rand();
  printf("a=\n");
  gf2_matrix_print(a);
  mm = gf2_build_matmul(a,0);
  if (mm) {
    printf("gf2_matmul build succeeded\n");
  } else return 0;

  b = gf2_matrix_alloc(16,16,field);
  for(i=0;i<16*16;i++) b->element[i]=rand();
  printf("b=\n");
  gf2_matrix_print(b);
  mm2 = gf2_build_mulmat(b,0);
  if (mm2) {
    printf("gf2_mulmat build succeeded\n");
  } else return 0;

  struct timeval tstart,tend;
  long long start,end;
  c=0;
  printf("normal mul:\n");
  gettimeofday(&tstart,0);
  for(i=0;i<TEST_ITER;i++)  c=gf2_matrix_mul(a,b,c);
  gettimeofday(&tend,0);
  if (c) {
    printf("Multiplied by gf2_matrix_mul\n");
    gf2_matrix_print(c);
    free(c);
    c=0;
  }
  start = tstart.tv_sec;
  start = start*1000000 + tstart.tv_usec;
  end = tend.tv_sec;
  end = end*1000000 + tend.tv_usec;
  printf("%lld microseconds\n",end-start);


  mm->accelerated_matmul = gf2_accel8_x86_64_matmul;
  mm->accelerated_matmul = gf2_accel4_x86_matmul;
  mm->accelerated_matmul = 0;
  mm->accelerated_matmul = gf2_accel16_sse_matmul;
  printf("mm->table = %llx\n",(long long unsigned)mm->table);
  c=0;
  printf("accelerated matmul:\n");
  gettimeofday(&tstart,0);
  for(i=0;i<TEST_ITER;i++) c = gf2_matrix_matmul(mm,b,c);
  gettimeofday(&tend,0);
  if (c) {
    printf("Multiplied by gf2_matrix_matmul\n");
    gf2_matrix_print(c);
    free(c);
    c=0;
  }
  start = tstart.tv_sec;
  start = start*1000000 + tstart.tv_usec;
  end = tend.tv_sec;
  end = end*1000000 + tend.tv_usec;
  printf("%lld microseconds\n",end-start);
  


  mm2->accelerated_mulmat = gf2_accel16_sse_mulmat;

  printf("accelerated mulmat:\n");
  gettimeofday(&tstart,0);
  for(i=0;i<TEST_ITER;i++) c = gf2_matrix_mulmat(a,mm2,c);
  gettimeofday(&tend,0);
  if (c) {
    printf("Multiplied by gf2_matrix_mulmat\n");
    gf2_matrix_print(c);
    free(c);
    c=0;
  }
  start = tstart.tv_sec;
  start = start*1000000 + tstart.tv_usec;
  end = tend.tv_sec;
  end = end*1000000 + tend.tv_usec;
  printf("%lld microseconds\n",end-start);

  free(a);
  free(b);

  int rows=1024;
  int bytes=128;
  
  a = gf2_matrix_alloc(4,bytes,field); /* 4 parity bytes for "bytes" data bytes */


  b = gf2_matrix_alloc(bytes,4,field); /* 4 parity bytes for "bytes" data bytes */

  for(i=0;i<bytes*4;i++) {
    a->element[i]=rand();
    b->element[i]=rand();

  }

  gf2_free_matmul(mm);
  gf2_free_mulmat(mm2);

  printf("a=\n");
  gf2_matrix_print(a);
  mm = gf2_build_matmul(a,0);
  if (mm) {
    printf("gf2_matmul build succeeded\n");
  } else return 0;

  printf("b=\n");
  gf2_matrix_print(b);
  mm2 = gf2_build_mulmat(b,0);
  if (mm2) {
    printf("gf2_mulmat build succeeded\n");
  } else return 0;



  d = gf2_matrix_alloc(bytes,rows,field);

  for(i=0;i<bytes*rows;i++) d->element[i]=rand();

  mm->accelerated_matmul = gf2_accel4_x86_matmul;
  printf("mm->table = %llx\n",(long long unsigned)mm->table);
  printf("accelerated matmul:\n");
  gettimeofday(&tstart,0);
  for(i=0;i<TEST_ITER;i++) c = gf2_matrix_matmul(mm,d,c);
  gettimeofday(&tend,0);
  if (c) {
    printf("Multiplied by gf2_matrix_matmul\n");
  }
  start = tstart.tv_sec;
  start = start*1000000 + tstart.tv_usec;
  end = tend.tv_sec;
  end = end*1000000 + tend.tv_usec;
  printf("%lld microseconds\n",end-start);
  printf("%d bytes total\n", bytes * rows * TEST_ITER);
  printf("%d result bytes\n", c->rows * c->columns * TEST_ITER);
  printf("%f data bytes/s\n", (float)(bytes * rows) * 1000000.0 * TEST_ITER / (float)(end-start)); 
  

  e = gf2_matrix_alloc(rows,bytes,field);
  for(i=0;i<bytes*rows;i++) e->element[i]=rand();


  mm2->accelerated_mulmat = gf2_accel4_x86_mulmat;
  printf("accelerated mulmat:\n");
  gettimeofday(&tstart,0);
  for(i=0;i<TEST_ITER;i++) h = gf2_matrix_mulmat(e,mm2,0);
  gettimeofday(&tend,0);
  if (h) {
    printf("Multiplied by gf2_matrix_mulmat\n");
  }
  start = tstart.tv_sec;
  start = start*1000000 + tstart.tv_usec;
  end = tend.tv_sec;
  end = end*1000000 + tend.tv_usec;
  printf("%lld microseconds\n",end-start);
  printf("%d bytes total\n", bytes * rows * TEST_ITER);
  printf("%d result bytes\n", h->rows * h->columns * TEST_ITER);
  printf("%f data bytes/s\n", (float)(bytes * rows) * 1000000.0 * TEST_ITER / (float)(end-start)); 

  f=gf2_matrix_mul(a,d,0);
  g=gf2_matrix_mul(e,b,0);

  for(i=0;i<rows*4;i++) {
    if (f->element[i] != c->element[i]) {
      printf("a*d != a*d at %i\n",i);
    }
    if (g->element[i] != h->element[i]) {
      printf("e*b != e*b at %i\n",i);
    }
  }


  printf("gf2 testing is complete.\n");
  return 0;
}
