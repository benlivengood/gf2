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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gf2.h"

/* Multiply two 7th degree polynomial elements of the field together
** using log and exp tables. The function can only take valid polynomials,
** so the only precondition is that field point to a valid gf2 structure that
** has been initialized, although in no case will the function crash
** if the structure has not been initialized */
uint8_t gf2_mul1(const gf2 *field,unsigned int a,unsigned int b)
{
  int unsigned i,j;
  if (a && b) { /* Both a and b must be nonzero to have valid logarithms */
    i=field->logtable[a]; /* Obtain logarithms of a and b, not that i and j */
    j=field->logtable[b]; /* range only from 0 to 255 */
    return field->exptable[i+j]; /* Return 2^(i+j), i+j is at most 510 */
  } else return 0; /* Return 0 if a or b was zero */
}

/* Multiply elements a and b using the 65536 byte multiplication table.
** Preconditions require that field points to a properly sized gf2 structure
** and that a and b are in the range 0-255. No error checking is performed */
uint8_t gf2_mul2(const gf2 *field,unsigned int a,unsigned int b)
{
  return field->multable[a+(b<<8)];
}

/* This function has no way to signal if b is 0, which is undefined. It
** returns a in this case, since the logtable has log(0)=255
** It calculates i=gf2_mulinv(field,a,b) such that i*b=a
** Preconditions: a and b in (0..255), field is a valid pointer to gf2 struct.
** Computes log(a/b) as log(a) - log(b), then exp(log(a)-log(b))
** Since the exp array is indexed from 0 to 511, log(a) must be adjusted to
** a value such that exp(log(a))=exp(log(a)+x), x>=255 so that log(b)<=x
** and the final array index is >=0. For exp(log(a))=exp(log(a)+x),
** exp(log(a))=exp(log(a)+x)) = exp(log(a))*exp(x), divide exp(log(a)) out to 
** get 1=exp(x) which means that x must be either 0 or 1 less than the order
** of the group by Fermat's theorem. In this case the group order is 256, so
** exp(log(a)-log(b) = exp(log(a)+255 - log(b))
*/
uint8_t gf2_mulinv(const gf2 *field,unsigned int a,unsigned int b)
{
  int i,j;
  if (a) { /* If a is zero the result of multiplication will be 0 for all b */
    i=field->logtable[a]+255; /* Add the order of our field to log(a) */
    j=field->logtable[b]; /* Logarithm of b; log(0)=0 in this table */
    return field->exptable[i-j]; /* Subtract logs to divide; exp(log(a/b)) */
  } else return 0; /* If a was zero, log tables can't handle, return it here */
}

/* Returns the inverse of a, without error checking 
** preconditions: a in (0..255), field a valid pointer to gf2 struct
** Computes inverse by exp(-log(a)) = 1 / exp(log(a))
** As above, the array index must be nonnegative, so compute 
** exp(-log(a))=exp(x-log(a)) for some positive x, and get 
** 1=exp(x) again, which equals 255 for the least nonnegative x.
*/
uint8_t gf2_inv(const gf2 *field,unsigned int a)
{
  return field->exptable[255 - field->logtable[a]];
}

/* Initalizes a gf2 structure for operations on a GF(2) field over an 7th
** degree polynomial defined by the low order 8 bits passed in poly.
** For bit b0..b7 (ordered LSB to MSB) the polynomial equals 
** x^8 + Sum(b_n * x^n). 
** Addition in the field is done with XOR. Doubling is done by left shift
** and XORing of the polynomial if the MSB was set.
** These two functions allow the multiplication tables to be built without
** having a function for general polynomial multiplication. First the
** logarithm and exponentiation tables are built, and then the large 256*256
** entry multiplication table can be generated.
** Preconditions: field is a valid pointer to a gf2 struct, poly represents
** an irreducible polynomial in GF(2). The irreducibility of the polynomial
** is tested straightforwardly by ensuring that every entry in the log
** table is unique. If the polynomial is not
** irreducible, 2^i will not generate the field of 256 possible polynomials,
** and therefore the logarithm table will have non-unique values. If the
** polynomial is irreducible, then 2^i will generate every element of the field
** for i=1 to 255. Once the algorithm ensures that all 255 nonzero values were
** generated, the polynomial is known to be irreducible and the log and exp
** tables are also complete. Actually, the loop progresses from 0 to 255, so
** there is one duplicate log entry written, namely log(x)=0 for both i=0 and
** i=255. This allows the exponential table to be built completely for all
** values from i=0 to 255, and the check for duplicate log values succeeds
** because the algorithm checks for a nonzero logarithm for 2^i. This
** test does not succeed for 2^0 or 2^255, because it is a zero value.
**
** Returns 0 on success, -1000+failing element on error
*/
int gf2_init(gf2 *field,char unsigned poly)
{
  int i,j;
  char unsigned element=1; /* Start with 2^0 = 1 */
  memset(field->logtable,0,256); /* Initialize the logtable for checking */
  field->logtable[0]=255; /* Invalid logarithm value will be 255, since
			  ** x^255 = x^0 over GF(2^8), so log(x)<255 */
  field->exptable[510]=1; /* x^(255*2) = 1, only used by mulinv in the
			  ** degenerate case of exp(1,0), which looks up
			  ** log(0)=255, adds 255 = 510, subtracts log(1)=0
			  ** and gets exptable[510] */
  for(i=0;i<255;i++) { /* i is the power to raise 2 to; element=2^i */
    field->exptable[i]=element; /* 2^i=element, so exp(i)=element */

    /* 2^i=exp(i), so exp(i+255)= exp(i)*exp(255) = exp(i)*2^255 = exp(i)*1
    ** by Fermat's theorem. */
    field->exptable[i+255]=element;

    /* Check to make sure this is a proper field. If any element is generated
    ** twice by 2^i, then the polynomial was reducible. Additionally,
    ** if 2^i=0 for any i, the polynomial was reducible.
    ** Since logtable was initialized to all zeros, if logtable[element]
    ** is nonzero the element was previously generated. Additionally, if
    ** element is zero then the polynomial is reducible, and logtable[0]=255 */
    if (field->logtable[element]!=0) return -1000+element;
    /* Check passed, log(2^i)=i */
    field->logtable[element]=i;
    /* If the MSB is set in the element, then 2*element is reduced by adding
    ** (XORing) element by the polynomial*2. The 8th degree is automatically
    ** reduced because a character only holds the low 8 coefficients, causing
    ** the x^8 term to drop off automatically. */

    if (element&0x80) element=(element<<1)^poly;
    else element<<=1;

    /* Fake code for generating tables based on 3^x instead of 2^x */
    /*    if (element&0x80) element^=(element<<1)^poly;
	  else element^=(element<<1);*/

  }
  /* Now that the log and exp tables are complete the multiplication table
  ** can be generated */

  memset(field->multable,0,256);
  for(i=1;i<256;i++) {
    int logi=field->logtable[i];
    uint8_t *exp=field->exptable+logi;
    for(j=1;j<i;j++) {
      int product=exp[field->logtable[j]];
      field->multable[i+(j<<8)]=product;
      field->multable[j+(i<<8)]=product;
    }
    field->multable[i<<8]=0;
    field->multable[i+(i<<8)]=field->exptable[logi+logi];
  }
  return 0; /* Successful creation of tables */
}


/* Exponentiation for x^y is performed as exp(log(x)*y) with the restriction
** that log(x)*y must fit in the exp table. This can be accomplished by
** finding a value z such that exp(z)=1 and dividing exp(log(x)*y) by
** exp(z) until log(x)*y is in (0..255). exp(255)=1 by Fermat's theorem, so
** exp(log(x)*y) = exp(log(x)*y)/n*exp(255) = exp(log(x)*y - n*255)
** = exp(log(x)*y mod 255)
*/
uint8_t gf2_exp(const gf2 *field,unsigned int x,unsigned int y)
{
  if (y==0) return 1; /* Anything to the power of zero is one */
  if (x) {
    x=field->logtable[x]; /* x=log(x) */
    return field->exptable[((x*y)%255)]; /* exp(log(x)*y mod 255) */
  } 
  return 0; /* If x is zero, any power is zero */
}

/* Preconditions: field is valid and a in [1..255]. Returns 0 if a=0 */
uint8_t gf2_log(const gf2 *field,unsigned int a)
{
  return field->logtable[a];
}

gf2_matrix *gf2_matrix_alloc(int unsigned rows,int unsigned columns,
			     const gf2 *field)
{
  gf2_matrix *tmp;
  int result;
  result = posix_memalign(&tmp,16,sizeof(gf2_matrix)+rows*columns);
  if (result!=0) return 0;
  tmp->rows=rows;
  tmp->columns=columns;
  tmp->field=field;
  return tmp;
}

gf2_matrix *gf2_matrix_add(const gf2_matrix *a, const gf2_matrix *b,
			   gf2_matrix *d)
{
  if (!a || !b || a->field != b->field ||
      a->rows != b->rows || a->columns != b->columns) return 0;
  if (d && ( d->field != a->field || d->rows != a->rows || 
	     d->columns != a->columns )) return 0;
  if (d==0) {
    d=gf2_matrix_alloc(a->rows, a->columns, a->field);
    if (d==0) return 0;
  }
  int i,s=a->rows*a->columns;
  for(i=0; i<s; i++) {
    d->element[i] = a->element[i] ^ b->element[i];
  }
  d->field=a->field;
  return d;
}

/* Calculate a*b by the matrix multiplication method.
** The field of a, b, and d (if nonzero) must be equal
** the number of columns in a must equal the number of rows in b.
** If d is a valid pointer to a matrix at with at least as many rows as a,
** and at least as many columns as b; or rather that the product of d's rows
** and columns is greater than or equal to the product of a's rows and b's
** columns.
** If d points to either a or b (or both) the function allocates temporary
** storage for the result and assigns the result to d after the computation
** is complete.
** If d is zero, a new gf2_matrix of appropriate dimensions is allocated and
** the result is placed in the new structure.
** The function returns the result, either d or the newly allocated matrix
** or returns 0 upon error.
*/
gf2_matrix *gf2_matrix_mul(const gf2_matrix *a,const gf2_matrix *b,
			   gf2_matrix *d)
{
  int i,j,k,di,ok,bi;
  uint8_t tmp;
  if (!a || !b || a->field != b->field ||
      a->columns != b->rows) return 0;
  if (d && ( d->field != a->field || d->rows*d->columns < a->rows*b->columns))
    return 0;
  if (d==0) { /* Allocate new matrix if d is null */
    d=gf2_matrix_alloc(a->rows, b->columns, a->field);
    if (d==0) return 0;
  }
  if (d==a || d==b) {
    uint8_t temp[a->rows*b->columns];
    for(i=0,di=0,bi=0; i<b->columns; i++,bi+=b->rows) { /* Loop on columns outer */
      for(j=0;j<a->rows;j++) { /* Loop on rows for the inner loop */
	tmp=0;
	for(k=0,ok=j; k<a->columns; k++,ok+=a->rows) {
	  tmp ^= GF2_MUL2(a->field->multable,a->element[ok],b->element[bi+k]);
	}
	temp[di++]=tmp;
      }
    }
    memcpy(d->element,temp,a->rows*b->columns);
  } else {
    for(i=0,di=0,bi=0; i<b->columns; i++,bi+=b->rows) { /* Loop on columns outer */
      for(j=0;j<a->rows;j++) { /* Loop on rows for the inner loop */
	tmp=0;
	for(k=0,ok=j; k<a->columns; k++,ok+=a->rows) {
	  tmp ^= GF2_MUL2(a->field->multable,a->element[ok],b->element[bi+k]);
	}
	d->element[di++]=tmp;
      }
    }
    d->rows=a->rows;
    d->columns=b->columns;
    d->field=a->field;
  }
  return d;
}



gf2_matrix *gf2_matrix_transpose(const gf2_matrix *a,gf2_matrix *d)
{
  if (!a) return 0;
  if (d && (d->field != a->field || d->rows*d->columns < a->rows*a->columns))
    return 0;
  if (d==0) { /* Allocate new matrix if d is null */
    d=gf2_matrix_alloc(a->columns, a->rows, a->field); /* transposed */
    if (d==0) return 0;
  }
  int i,j,ao,desto;
  if (d==a) {
    uint8_t temp[a->columns*a->rows];
    int unsigned rt;
    for(i=0,ao=0; i<a->columns; i++,ao += a->rows) {
      for(j=0,desto=0; j<a->rows; j++,desto += a->columns) {
	temp[desto+i]=a->element[ao+j];
      }
    }
    memcpy(d->element,temp,a->columns*a->rows);
    rt=a->rows;
    d->rows=a->columns; /* Same matrix, so swap columns and rows */
    d->columns=rt;
  } else {
    d->rows=a->columns;
    d->columns=a->rows;
    d->field=a->field;
    for(i=0,ao=0; i<a->columns; i++,ao += a->rows) {
      for(j=0,desto=0; j<a->rows; j++,desto += a->columns) {
	d->element[desto+i]=a->element[ao+j];
      }
    }
  }
  return d;
}

gf2_matrix *gf2_gaussian_elimination(const gf2_matrix *a,gf2_matrix *d)
{
  gf2_matrix *matrix;
  if (!a) return 0;

  matrix=gf2_matrix_alloc(a->rows,a->columns,a->field);
  if (matrix==0) return 0;


  int rows=a->rows,cols=a->columns;
  /* If this matrix has more columns than rows, transpose it to solve it */
  if (a->rows<a->columns) {
    cols=a->rows;
    rows=a->columns;
    if (!gf2_matrix_transpose(a,matrix)) {
      free(matrix);
      return 0;
    }
  } else memcpy(matrix->element,a->element,rows*cols);



  /* First pass of gaussian elimination. Make left column equal to 1 x x ...
  ** and then make the second column equal to 0 1 x ... and continue until
  ** the last row is ... 0 0 1
  ** Original code worked row by row. With change to column major matrices
  ** the algorithm works best column by column
  */

  int i,j,k,oj,ok;
  char unsigned cur;
  const char unsigned *mt=a->field->multable;
  for(k=0,ok=0;k<cols;k++,ok+=rows) { /* For every column */
    for(j=0,oj=0;j<k;j++,oj+=rows) { /* For every entry before the diagonal */
      cur=matrix->element[ok+j]; /* Multiplicand entry we're using */
      for(i=0;i<rows;i++) { /* For every entry in the column */
	matrix->element[ok+i]^=GF2_MUL2(mt,cur,matrix->element[oj+i]); /* subtract other row */
      }      
    }
    if (matrix->element[ok+k]==0) {
      free(matrix);
      return 0; /* This entry is invalid, matrix is non-invertible */
    }
    cur=gf2_inv(a->field,matrix->element[ok+k]); /* Diagonal entry from matrix*/
    for(i=k;i<rows;i++) { /* For every entry in the row (not before diag) */
      matrix->element[ok+i]=GF2_MUL2(mt,cur,matrix->element[ok+i]); /* Multiply by constant */
    }
  }
  /* Second pass of gaussian elimination from the bottom up */
  for(k=cols-2,ok=(cols-2)*rows;k>=0;k--,ok-=rows) {
    for(j=k+1,oj=ok+rows;j<cols;j++,oj+=rows) {
      cur=matrix->element[ok+j]; /* Multiplicand entry we're using */
      for(i=0;i<rows;i++) { /* For every entry in the row */
	matrix->element[ok+i]^=GF2_MUL2(mt,cur,matrix->element[oj+i]); /* subtract other row */
      }      
    }
  }

  /* If the original matrix had more columns than rows, transpose the
  ** solution into the result destination */
  if (d) {
    if (a->rows<a->columns) {
      gf2_matrix_transpose(matrix,d);
    } else memcpy(d->element,matrix->element,rows*cols);
    d->field=a->field;
    free(matrix);
    return d;
  } else {
    if (a->rows<a->columns) {
      gf2_matrix_transpose(matrix,matrix);
    }
    return matrix;
  }
}

void gf2_matrix_print(const gf2_matrix *a)
{
  int i,j,io;
  if (!a) return;
  for(j=0; j<a->rows; j++) {
    for(i=0,io=0; i<a->columns; i++,io+=a->rows) {
      printf("%.2x ",a->element[io+j]);
    }
    printf("\n");
  }
}

gf2_matrix *gf2_matrix_set_zero(gf2_matrix *a)
{
  if (!a) return 0;
  memset(a->element,0,a->rows*a->columns);
  return a;
}

gf2_matrix *gf2_matrix_set_identity(gf2_matrix *a)
{
  int i;
  if (!a) return 0;
  memset(a->element,0,a->rows*a->columns);
  for(i=0;i<a->rows && i<a->columns;i++) {
    a->element[i+i*a->rows]=0x01;
  }
  return a;
}

gf2_matrix *gf2_matrix_copy(const gf2_matrix *a,gf2_matrix *d,
			    int unsigned fromrow,int unsigned fromcol,
			    int unsigned torow,int unsigned tocol,
			    int unsigned rows,int unsigned cols)
{
  int i,j,fromj,toj;
  if (!a) return 0;
  if (fromrow<0 || fromcol<0 || torow<0 || tocol<0) return 0;
  if (fromrow+rows>a->rows || fromcol+cols>a->columns) return 0;
  if (d && (a->field != d->field || 
	    d->rows*d->columns < (torow+rows)*(tocol+cols)))
    return 0;
  if (!d) {
    d=gf2_matrix_alloc(torow+rows,tocol+cols,a->field);
    if (!d) return 0;
  } else {
    if (torow+rows>d->rows || tocol+cols>d->columns) return 0;
  }
  if (d==a) {
    uint8_t temp[d->rows*d->columns];
    memcpy(temp,d->element,d->rows*d->columns);
    for(j=0,fromj=fromcol*a->rows,toj=tocol*d->rows ; j<cols ;
	j++,fromj+=a->rows,toj+=d->rows) {
      for(i=0;i<rows;i++) {
	temp[i+torow + toj] = a->element[i+fromrow + fromj];
      }
    }
    memcpy(d->element,temp,d->rows*d->columns);
  } else {
    for(j=0,fromj=fromcol*a->rows,toj=tocol*d->rows ; j<cols ;
	j++,fromj+=a->rows,toj+=d->rows) {
      for(i=0;i<rows;i++) {
	d->element[i+torow + toj] = a->element[i+fromrow + fromj];
      }
    }
    d->field=a->field;
  }
  return d;
}


gf2_matrix *gf2_matrix_mul_scalar(const gf2_matrix *a,uint8_t scalar,
				  gf2_matrix *d)
{
  int i;
  const uint8_t *multable;
  if (!a) return 0;
  if (!d) {
    d=gf2_matrix_alloc(a->rows,a->columns,a->field);
    if (!d) return 0;
  }
  multable=a->field->multable+(((int unsigned)scalar)<<8);
  for(i=0;i<a->rows*a->columns;i++) {
    d->element[i]=multable[a->element[i]];
  }
  d->field=a->field;
  return d;
}


gf2_matmul *gf2_build_matmul(const gf2_matrix *src,gf2_matmul *dest)
{
  unsigned int table_rows,byte_size;
  uint8_t *tmptable;
  if (src->rows <= 4) table_rows=4; /* 4 byte table entries */
  else if (src->rows <= 8) table_rows=8; /* or 8 byte entries */
  else table_rows = (src->rows + 0xf) & ~(0xf); /* or a multiple of 16 */

  byte_size=256 * src->columns * table_rows; /* Number of bytes in the table */
  
  if (dest) {
    if (dest->table) {
      if (dest->columns * dest->table_rows * 256 < byte_size) {
	int result;
	result = posix_memalign(&tmptable,16,table_rows * 256 * src->columns);
	if (result !=0) return 0;
	free(dest->table);
	dest->table = tmptable;
      }
    }
  } else {
    int result;


    dest = (gf2_matmul *)malloc(sizeof(gf2_matmul));
    if (dest == 0) return 0;
    /*    printf("alloc: %d\n",table_rows * 256 * src->columns);*/
    result = posix_memalign(&(dest->table),16,table_rows * 256 * src->columns);
    if (result != 0) {
      free(dest);
      return 0;
    }
  }
  dest->columns = src->columns;
  dest->rows = src->rows;
  dest->table_rows = table_rows;
  dest->field = src->field;
  dest->accelerated_matmul = 0; /* By default there is no accelerated function */

  /*  printf("Source for matmul: (%d,%d)\n",src->rows,src->columns);
      printf("Dest for matmul: (%d,%d)\n",dest->rows,dest->columns);*/

  /*  printf("columns=%d  rows=%d  table_rows=%d  field=%x\n",dest->columns,
      dest->rows,dest->table_rows,dest->field);*/

  int i,j,k,srccol;
  const uint8_t *mul;
  uint8_t *tp = dest->table;
  for(i=0;i<src->columns;i++) { /* For each column in the matrix */
    srccol = i * src->rows;
    /*    printf("i=%d  srccol=%d\n",i,srccol); */
    for(j=0;j<256;j++) { /* For each of the 256 possible GF(2) elements */
      /*      printf("j=%d\n",j); */
      mul = dest->field->multable + j*256; /* Load the multiplication table */
      for(k=0;k<dest->rows;k++) *tp++ = mul[src->element[srccol + k]]; /* multiply column */
      for(;k<table_rows;k++) *tp++ = 0; /* Zero fill any remaining bytes */
    }
  }
  return dest;
}

/* This is a generic function that does bytewise table application
** The destination matrix will have src->rows columns and mm->columns rows.
** dest CAN NOT point to the same matrix as src
*/

gf2_matrix *gf2_matrix_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			      gf2_matrix *dest)
{
  if (mm->columns != src->rows) return 0;
  
  if (dest) {
    if (dest->rows * dest->columns < mm->table_rows * src->columns) return 0;
    dest->rows = mm->table_rows;
    dest->columns = src->columns;
  } else {
    dest = gf2_matrix_alloc(mm->table_rows,src->columns,
			    src->field);
    if (dest == 0 ) return 0;
  }
  
  /*  printf("src size: (%d,%d)\n",src->rows,src->columns);*/
  
  if (mm->accelerated_matmul) mm->accelerated_matmul(mm,src,dest);
  else {
    
    int i,j,k,oj,uint32rows,oi;
    uint32_t *table,*dest32;
    uint32rows = mm->table_rows >> 2;
    
    dest32 = (uint32_t *)dest->element;
    for(j=0,oj=0;j<src->columns;j++,oj+=src->rows) { /* Outer loop over the rows */
      table=(uint32_t *)(mm->table) + src->element[oj]*uint32rows;
      for(k=0;k<uint32rows;k++) {
	dest32[k] = table[k];
      }
      for(oi=256*uint32rows,i=1;i<src->rows;i++,oi += 256*uint32rows) {
	table = (uint32_t *)(mm->table) + oi + src->element[oj+i]*uint32rows;
	for(k=0;k<uint32rows;k++) {
	  dest32[k] ^= table[k];
	}
      }
      dest32 += uint32rows;
    }
  }
  return dest;
}


void gf2_accel4_x86_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			      gf2_matrix *dest)
{
  int i,j,oj;
  uint32_t *table,*dest32,tmp;
  
  dest32 = (uint32_t *)dest->element;
  for(j=0,oj=0;j<src->columns;j++,oj+=src->rows) { /* Outer loop over the rows */
    table=(uint32_t *)(mm->table);
    tmp = table[src->element[oj]];
    for(i=1;i<src->rows;i++) {
      table += 256;
      tmp ^= table[src->element[oj+i]];
    }
    *dest32++ = tmp;
  }
  return;
}


void gf2_accel8_x86_64_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			      gf2_matrix *dest)
{
  int i,j,oj;
  uint64_t *table,tmp;
   
  for(j=0,oj=0;j<src->columns;j++,oj+=src->rows) { /* Outer loop over the rows */
    table=(uint64_t *)(mm->table);
    tmp = table[src->element[oj]];
    for(i=1;i<src->rows;i++) {
      table += 256;
      tmp ^= table[src->element[oj + i]];
    }
    ((uint64_t *)dest->element)[j] = tmp;
  }
  
  return;
}


void gf2_accel16_sse_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			      gf2_matrix *dest)
{
  int i,j,oj;
  v4si *table,tmp;
   
  for(j=0,oj=0;j<src->columns;j++,oj+=src->rows) { /* Outer loop over the rows */
    table=(v4si *)(mm->table);
    tmp = table[src->element[oj]];
    for(i=1;i<src->rows;i++) {
      table += 256;
      tmp ^= table[src->element[oj + i]];
    }
    ((v4si *)dest->element)[j] = tmp;
  }
  
  return;
}


void gf2_free_matmul(gf2_matmul *mm)
{
  if (mm->table) free(mm->table);
  free(mm);
}

/* mulmat:
** Multiply a normal matrix by an accelerated matrix
**
**  a00 a01 a02 a03 a04             b00 b01     m00 m01
**  a10 a11 a12 a13 a14   *         b10 b11     m10 m11
**  a20 a21 a22 a23 a24             b20 b21  =  m20 m21
**  a30 a31 a32 a33 a34             b30 b31     m30 m31
**  a40 a41 a42 a43 a44             b40 b41     m40 m41
**  a50 a51 a52 a53 a54                         m50 m51
**
**
** but the matrix is built as:
**
**   m00 m10 m20 m30 m40 m50
**   m01 m11 m21 m31 m41 m51
**
** because that's how the acceleration structure works. We need to
** transpose the result at the end of the multiplication step
**
*/

/* Start of mulmat functions */
gf2_mulmat *gf2_build_mulmat(const gf2_matrix *src,gf2_mulmat *dest)
{
  unsigned int table_columns,byte_size;
  uint8_t *tmptable;
  if (src->columns <= 4) table_columns=4; /* 4 byte table entries */
  else if (src->columns <= 8) table_columns=8; /* or 8 byte entries */
  else table_columns = (src->columns + 0xf) & ~(0xf); /* or a multiple of 16 */

  byte_size=256 * src->rows * table_columns; /* Number of bytes in the table */
  
  if (dest) {
    if (dest->table) {
      if (dest->rows * dest->table_columns * 256 < byte_size) {
	int result;
	result = posix_memalign(&tmptable,16,table_columns * 256 * src->rows);
	if (result !=0) return 0;
	free(dest->table);
	dest->table = tmptable;
      }
    }
  } else {
    int result;
    
    
    dest = (gf2_mulmat *)malloc(sizeof(gf2_mulmat));
    if (dest == 0) return 0;
    /*    printf("alloc: %d\n",table_rows * 256 * src->columns);*/
    result = posix_memalign(&(dest->table),16,table_columns * 256 * src->rows);
    if (result != 0) {
      free(dest);
      return 0;
    }
  }
  dest->columns = src->columns;
  dest->rows = src->rows;
  dest->table_columns = table_columns;
  dest->field = src->field;
  dest->accelerated_mulmat = 0; /* By default there is no accelerated function */

  /*  printf("Source for matmul: (%d,%d)\n",src->rows,src->columns);
      printf("Dest for matmul: (%d,%d)\n",dest->rows,dest->columns);*/

  /*  printf("columns=%d  rows=%d  table_rows=%d  field=%x\n",dest->columns,
      dest->rows,dest->table_rows,dest->field);*/

  int i,j,k,srccol;
  const uint8_t *mul;
  uint8_t *tp = dest->table;
  for(i=0;i<src->rows;i++) { /* For each row in the matrix */
    /*    printf("i=%d  srccol=%d\n",i,srccol); */
    for(j=0;j<256;j++) { /* For each of the 256 possible GF(2) elements */
      /*      printf("j=%d\n",j); */
      mul = dest->field->multable + j*256; /* Load the multiplication table */
      for(srccol=i,k=0;k<dest->columns;k++,srccol+=src->rows) *tp++ = mul[src->element[srccol]]; /* multiply column */
      for(;k<table_columns;k++) *tp++ = 0; /* Zero fill any remaining bytes */
    }
  }
  return dest;
}

/* This is a generic function that does bytewise table application
** The destination matrix will have src->columns columns and mm->rows rows.
** dest CAN NOT point to the same matrix as src
*/

gf2_matrix *gf2_matrix_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			      gf2_matrix *dest)
{
  if (mm->rows != src->columns) return 0;
  
  if (dest) {
    if (dest->rows * dest->columns < mm->table_columns * src->rows) return 0;
    /*    dest->rows = src->rows;
	  dest->columns = mm->table_columns;*/
    /* We're building the transpose of the real matrix at this point */
    dest->rows = mm->table_columns;
    dest->columns = src->rows;
  } else {
    /*    dest = gf2_matrix_alloc(src->rows,mm->table_columns,
	  src->field);*/
    /* We're building the transpose of the real matrix at this point */
    dest = gf2_matrix_alloc(mm->table_columns,src->rows,
			    src->field);
    if (dest == 0 ) return 0;
  }
  
  /*  printf("src size: (%d,%d)\n",src->rows,src->columns);*/
  
  if (mm->accelerated_mulmat) mm->accelerated_mulmat(src,mm,dest);
  else {
    /* Multiplying a matrix by an accelerated matrix means that going across
    ** the row in the source matrix pulls bytes from non-adjacent addresses.
    ** 
    ** The solution is to take n rows at a time, and n columns from the
    ** destination matrix and work on them all at the same time.
    ** For the default 32 bit function we'll do 4 rows/columns at a time.
    */
    int i,j,k,oj,uint32columns = mm->table_columns >> 2;
    uint32_t *table,*dest32;
    dest32 =(uint32_t *)dest->element;
    for(i=0;i<src->rows;i++) { /* Work 1 row at a time */
      table = (uint32_t *)(mm->table) + src->element[i]*uint32columns;
      for(k=0;k<uint32columns;k++) {
	dest32[k]= table[k];
      }
      for(oj=256*uint32columns,j=1;j<src->columns;j++,oj+=256*uint32columns) {
	table = (uint32_t *)(mm->table) + oj + src->element[i+j*src->rows]*uint32columns;
	for(k=0;k<uint32columns;k++) {
	  dest32[k] ^= table[k];
	}
      }
      dest32 += uint32columns;
    }
    gf2_matrix_transpose(dest,dest);
  }
  return dest;
}


void gf2_accel4_x86_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			      gf2_matrix *dest)
{
  
  int di,i,j,oj;
  int off1,off2,off3;

  uint32_t *table,tmp;

  uint32_t tmp1,tmp2,tmp3,w0,w1,w2,w3;

  uint32_t *dest32;
  dest32 =(uint32_t *)dest->element;
    
  dest->columns=mm->columns;
  dest->rows=src->rows;

  off1=dest->rows>>2;
  off2=off1*2;
  off3=off1*3;
  
  for(di=0,i=0;i<src->rows;i+=4,di++) { /* Work 4 rows at a time */
    table = (uint32_t *)(mm->table);
    tmp = table[src->element[i]];
    tmp1 = table[src->element[i+1]];
    tmp2 = table[src->element[i+2]];
    tmp3 = table[src->element[i+3]];
    
    for(j=1,oj=src->rows;j<src->columns;j++,oj+=src->rows) {
      table+=256;
      tmp ^= table[src->element[i+oj]];
      tmp1 ^= table[src->element[i+oj+1]];
      tmp2 ^= table[src->element[i+oj+2]];
      tmp3 ^= table[src->element[i+oj+3]];
    }

    /* transpose tmp .. tmp3 into w0 .. w3
    **
    ** tmp  tmp1 tmp2 tmp3
    ** w0_0 w0_1 w0_2 w0_3
    ** w1_0 w1_1 w1_2 w1_3 
    ** w2_0 w2_1 w2_2 w2_3
    ** w3_0 w3_1 w3_2 w3_3  where w<index>_<byte> and byte from 0 to 3
    **
    */
    w0 = (tmp&0xff) | ((tmp1&0xff)<<8) | ((tmp2&0xff)<<16) | (tmp3<<24);
    w1 = ((tmp>>8)&0xff) | (tmp1&0xff00) | ((tmp2&0xff00)<<8) | ((tmp3&0xff00)<<16);
    w2 = ((tmp>>16)&0xff) | ((tmp1>>8)&0xff00) | (tmp2&0xff0000) | ((tmp3&0xff0000)<<8);
    w3 = (tmp>>24) | ((tmp1&0xff000000)>>16) | ((tmp2&0xff000000)>>8) | (tmp3&0xff000000);
    
    dest32[di]=w0;
    dest32[di+off1]=w1;
    dest32[di+off2]=w2;
    dest32[di+off3]=w3;
  }
}







void gf2_accel8_x86_64_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			      gf2_matrix *dest)
{
}


/*
** After the table lookups the resulting sse registers will be transposed from the
** result we want to store. It would be nice to get the results we want into memory in
** the right order using SSE operations 
**
** if we calculate:
**
**   xmm0 xmm1 xmm2 xmm3
**   p00  p10  p20  p30
**   p01  p11  p21  p31
**   p02  p12  p22  p32
**   p03  p13  p23  p33
**   p04  p14  p24  p34
**   p05  p15  p25  p35
**   p06  p16  p26  p36
**   p07  p17  p27  p37
**   p08  p18  p28  p38
**   p09  p19  p29  p39
**   p0a  p1a  p2a  p3a
**   p0b  p1b  p2b  p3b
**   p0c  p1c  p2c  p3c
**   p0d  p1d  p2d  p3d
**   p0e  p1e  p2e  p3e
**   p0f  p1f  p2f  p3f
**
** We need to put it in the form of
**
**   p00  p01  p02  p03  p04  p05  p06  p07  p08  p09  p0a  p0b  p0c  p0d  p0e  p0f
**   p10  p11  p12  p13  p14  p15  p16  p17  p18  p19  p1a  p1b  p1c  p1d  p1e  p1f
**   p20  p21  p22  p23  p24  p25  p26  p27  p28  p29  p2a  p2b  p2c  p2d  p2e  p2f
**   p30  p31  p32  p33  p34  p35  p36  p37  p38  p39  p3a  p3b  p3c  p3d  p3e  p3f
**
**
**   d00 d01 d02 d03     m00 m01 m02 m03   p00 p01 p02 p03
**   d10 d11 d12 d13     m10 m11 m12 m13   p10 p11 p12 p13
**   d20 d21 d22 d23     m20 m21 m22 m23   p20 p21 p22 p23
**   d30 d31 d32 d33  *  m30 m32 m32 m33 = p30 p31 p32 p33
**   d40 d41 d42 d43                       p40 p41 p42 p43
**   d50 d51 d52 d53                       p50 p51 p52 p53
**   d60 d61 d62 d63                       p60 p61 p62 p63
**   d70 d71 d72 d73                       p70 p71 p72 p73
**   d80 d81 d82 d83                       p80 p81 p82 p83
**   d90 d91 d92 d93                       p90 p91 p92 p93
**   da0 da1 da2 da3                       pa0 pa1 pa2 pa3
**   db0 db1 db2 db3                       pb0 pb1 pb2 pb3
**   dc0 dc1 dc2 dc3                       pc0 pc1 pc2 pc3
**   dd0 dd1 dd2 dd3                       pd0 pd1 pd2 pd3
**   de0 de1 de2 de3                       pe0 pe1 pe2 pe3
**   df0 df1 df2 df3                       pf0 pf1 pf2 pf3
**
**
**   p00 = (d00*m00 + d01*m10 + d02*m20 + d03*m30)
**   p10 = (d10*m00 + d11*m10 + d12*m20 + d13*m30)
**   p20 = (d20*m00 + d21*m10 + d22*m20 + d23*m30)
**
**
**
*/

void gf2_accel16_sse_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			      gf2_matrix *dest)
  {
  int i,j,off;
  v4si *table,*dest128,tmp0,tmp1,tmp2,tmp3;
  dest128 =(v4si *)dest->element;
  for(i=0;i<src->rows;i+=4) { /* Work 4 rows at a time */
    table = (v4si *)mm->table;
    tmp0 = table[src->element[i]];
    tmp1 = table[src->element[i+1]];
    tmp2 = table[src->element[i+2]];
    tmp3 = table[src->element[i+3]];
    table+=256;
    for(j=1,off=i+src->rows;j<src->columns;j++,off+=src->rows,table+=256) {
      tmp0 ^= table[src->element[off]];
      tmp1 ^= table[src->element[off+1]];
      tmp2 ^= table[src->element[off+2]];
      tmp3 ^= table[src->element[off+3]];
    }
    *dest128++ = tmp0;
    *dest128++ = tmp1;
    *dest128++ = tmp2;
    *dest128++ = tmp3;
  }
  gf2_matrix_transpose(dest,dest);
  return;

}


void gf2_free_mulmat(gf2_mulmat *mm)
{
  if (mm->table) free(mm->table);
  free(mm);
}









gf2_matrix *gf2_matrix_orthonormalize(const gf2_matrix *a,
				      gf2_matrix *d)
{




  return 0;
}
