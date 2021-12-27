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

#ifndef GF2_H
#define GF2_H

#include <stdint.h>

/* SSE 16-byte data type */
typedef int v4si __attribute__ ((vector_size(16)));


/* Slower multiplicatoin function requiring less memory */
#define GF2_MUL1(exp,log,a,b) (((a)&&(b))?((exp)[(log)[(a)]+(log)[(b)]]):0)

/* Faster multiplication function relying on large (64K) lookup table */
#define GF2_MUL2(mul,a,b) ((mul)[(a)+((b)<<8)])

/* Logarithm base 0x02 over GF(2^8) from table lookup */
#define GF2_LOG(field,a) ((field)->logtable[a])

/* Exponentiation of 0x02 over GF(2^8) from table lookup */
#define GF2_EXP(field,a) ((field)->exptable[a])

/* This structure contains the bytewise multiplication table lookup table.
** Multable contains i*j for multable[i*256+j]
** The exponential table contains 2^n for n=0 to 511
** The logtable which contains log_2(n) for 0 through 255, with the undefined
** value log_2(0) in the first element of the array and equal to zero.
** The polynomial is a bitwise representation of the lowest 8 bits of an
** implicit 8th degree GF(2) polynomial. Upon initialization the polynomial
** is tested for primality and thus for a true field over all 256 values.
** No destructor is needed for the structure because all values are stored
** locally and not allocated from the heap.
*/

typedef struct gf2 {
  uint8_t multable[65536];
  uint8_t logtable[256];
  uint8_t exptable[511]; /* log(x)<255, so 254*2=510 entries maximum */
  uint8_t polynomial;
} gf2;

/* Column major matrix structure. columns*rows+sizeof(int unsigned)*2+
** sizeof(const gf2 *) bytes, with a11 being the first byte, a21 being
** the second byte, etc.
*/

typedef struct gf2_matrix {
  int unsigned columns;
  int unsigned rows;
  const gf2 *field;
  uint8_t element[];
} gf2_matrix;

/* gf2_matmul is an acceleration structure that can be used to speed up
** matrix multiplications. It consists of an array of premultiplied columns of
** the matrix for each column position.
** Each premultiplied column has 256*rows entries, one entry of rows size for
** each of the elements of the GF2 field multiplied by the matrix column it
** represents.

 (m00 m01 m02 m03 m04 m05)   ( d0 )
 (m10 m11 m12 m13 m14 m15)   ( d1 )
                             ( d2 )  =  ( p0 )
                             ( d3 )     ( p1 )
                             ( d4 )
                             ( d5 )

  For instance, the above matrix would have 6 premultiplied tables, each with
256 2-byte entries. Then when multiplying M by D the result is calculated from
the GF(2) sum of the table entries for d0 through d5.

  The actual number of rows in the table is equal to the maximum of the
rows in the matrix and the least power of 2 greater than or equal to the number
of rows. This keeps the tables aligned, and any unused bytes in the table are
set to zero, this maximum is called table_rows in the structure.

The organization of the table is [columns][elements][table_rows], so that
all bytes in a row are contiguous, and all rows are continguous within a column
	
It is relatively important that table be aligned if SSE instructions are used.
		
accelerated_mul is a function pointer that can be called if it's nonzero.
It does FAST multiplication with the acceleration tables, using MMX and SSE
if possible.
     
**
*/

typedef struct gf2_matmul {
  int unsigned columns;
  int unsigned rows;
  int unsigned table_rows;
  const gf2 *field;
  uint8_t *table;
  void (*accelerated_matmul)(const struct gf2_matmul *mm,const gf2_matrix *src,
			  gf2_matrix *dest);
} gf2_matmul;

/* Same type of structure as above, but this one is for multiplying a normal
** matrix by an accelerated matrix, which makes all the setup and operations
** essentially the transpose of gf2_matmul.
*/
typedef struct gf2_mulmat {
  int unsigned columns;
  int unsigned rows;
  int unsigned table_columns;
  const gf2 *field;
  uint8_t *table;
  void (*accelerated_mulmat)(const gf2_matrix *src,const struct gf2_mulmat *mm,
			  gf2_matrix *dest);
} gf2_mulmat;


/* Initialize a gf2 field structure given the lowest 8 coefficients of an 
** implicit 8th degree polynomial over GF(2). The polynomial is validated 
** to produce a field over all 256 elements by a simple test
*/
int gf2_init(gf2 *field,uint8_t poly);

/* Multiply two elements of GF(2) over our polynomial using log and exp
** tables. */
unsigned char gf2_mul1(const gf2 *field,unsigned int a,unsigned int b);

/* Multiply two elements of GF(2) over our polynomial using big multiplication
** table. This is faster on newer processors */
unsigned char gf2_mul2(const gf2 *field,unsigned int a,unsigned int b);

/* gfs_inv returns the inverse of a, unless a is zero in which case zero
** is also returned. Should this return an integer with -1 for undefined? */
unsigned char gf2_inv(const gf2 *field,unsigned int a);

/* gf2_mulinv multiplies a by the inverse of b, equiv. to a/b. b=0 == error */
unsigned char gf2_mulinv(const gf2 *field,unsigned int a,unsigned int b);

/* Returns x^y in the polynomial field. x is a field element, y an integer */
unsigned char gf2_exp(const gf2 *field,unsigned int x,unsigned int y);

/* Returns log(a) in the polynomial field */
unsigned char gf2_log(const gf2 *field,unsigned int a);

/* Allocate a rows by columns matrix using malloc, and pointing to field */
gf2_matrix *gf2_matrix_alloc(int unsigned rows,int unsigned columns,
			     const gf2 *field);

/* Calculate a+b by pairwise sum of the elements of a and b, e.g. a11 + b11,
** a23 + b23, etc.
** The field of a, b, and d (if nonzero) must be equal
** a and b must have the same dimensions.
** If d is a valid pointer to a matrix of the same dimensions as a and b
** the sum is placed in d.
** If d points to either a or b (or both) the function still behaves correctly.
** If d is zero, a new gf2_matrix is allocated and the sum placed in the new
** matrix.
** The function returns the destination matrix, either d or the newly
** allocated matrix, or returns 0 upon error.
** Addition is over GF(2) (XOR).
** 
*/
gf2_matrix *gf2_matrix_add(const gf2_matrix *a,const gf2_matrix *b,
			   gf2_matrix *d);

/* Calculate a*b by the matrix multiplication method.
** The field of a, b, and d (if nonzero) must be equal
** the number of columns in a must equal the number of rows in b.
** If d is a valid pointer to a matrix with the same number of rows as
** a and the same number of columns as b, then the result is placed in d.
** If d points to either a or b (or both) the function allocates temporary
** storage for the result and assigns the result to d after the computation
** is complete.
** If d is zero, a new gf2_matrix of appropriate dimensions is allocated and
** the result is placed in the new structure.
** The function returns the result, either d or the newly allocated matrix
** or returns 0 upon error.
*/
gf2_matrix *gf2_matrix_mul(const gf2_matrix *a,const gf2_matrix *b,
			   gf2_matrix *d);




/* Performs gaussian elimination on the matrix a, storing the result in
** dest. Both matricies must be rows by columns matrices. If rows>columns,
** the upper square of the matrix is solved by elementary column operations,
** otherwise row operations are performed.
** Returns 0 if gaussian elimination can be performed.
** Returns 1 for general errors. (rows, columns, a, or dest invalid)
** Returns 2 if the upper left square matrix is not linearly independent,
** which is detected by the inability to find an appropriate row to
** use in the next step of gaussian elimination.
*/
gf2_matrix *gf2_gaussian_elimination(const gf2_matrix *a,gf2_matrix *d);

/* Transpose the matrix a and return the result. If d is null, a new
** matrix is allocated and the result stored in it, otherwise the result
** is stored in d. Either d or the new matrix is returned, or zero on error.
** The field of a and d (if nonzero) must be equal
** the number of columns in a must equal the number of rows in d, and vice
** versa.
** If d is a valid pointer to a matrix, then the result is placed in d.
** If d points to a the function allocates temporary
** storage for the result and assigns the result to d after the computation
** is complete.
** If d is zero, a new gf2_matrix of appropriate dimensions is allocated and
** the result is placed in the new structure.
** The function returns the result, either d or the newly allocated matrix
** or returns 0 upon error.
*/

gf2_matrix *gf2_matrix_transpose(const gf2_matrix *a,gf2_matrix *d);
/*void gf2_matrix_transpose(const uint8_t *a,
			  uint8_t *dest,int unsigned arows,
			  int unsigned acols);*/

/* Print the matrix */
void gf2_matrix_print(const gf2_matrix *a);

/* Zero out the matrix */
gf2_matrix *gf2_matrix_set_zero(gf2_matrix *a);

/* Set the matrix to the identity matrix, or as close as possible for
** non-square matrices */
gf2_matrix *gf2_matrix_set_identity(gf2_matrix *a);

/* Copy a part of one matrix to a location in another matrix */
gf2_matrix *gf2_matrix_copy(const gf2_matrix *a,gf2_matrix *d,
			    int unsigned fromrow,int unsigned fromcol,
			    int unsigned torow,int unsigned tocol,
			    int unsigned rows,int unsigned cols);

/* Multiply a matrix by a scalar. Each element is multiplied by the 
** scalar value.
** If d is null, a new matrix is allocated and returned.
** null is returned on error
*/
gf2_matrix *gf2_matrix_mul_scalar(const gf2_matrix *a,uint8_t scalar,
				  gf2_matrix *d);


/* Build a gf2_matmul struct from an existing matrix. If dest is non-zero and
** has enough space (as calculated by 256*table_columns*rows), it is filled
** in with the new tables, otherwise zero is returned. If dest is zero,
** a new gf2_matmul structure is malloced. If that fails, zero is returned.
** If everything succeeds, a pointer to gf2_matmul is returned.
 */
gf2_matmul *gf2_build_matmul(const gf2_matrix *src,gf2_matmul *dest);

void gf2_accel8_x86_64_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			 gf2_matrix *dest);
void gf2_accel4_x86_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			   gf2_matrix *dest);
void gf2_accel16_sse_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			    gf2_matrix *dest);

/* Calculate mm*src from a gf2_matmul acceleration structure for mm.
** Destination CAN NOT be the same as source.
*/
gf2_matrix *gf2_matrix_matmul(const gf2_matmul *mm,const gf2_matrix *src,
			      gf2_matrix *dest);

/* Free a gf2_matmul structure */
void gf2_free_matmul(gf2_matmul *mm);

gf2_mulmat *gf2_build_mulmat(const gf2_matrix *src,gf2_mulmat *dest);

void gf2_accel8_x86_64_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			 gf2_matrix *dest);
void gf2_accel4_x86_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			   gf2_matrix *dest);
void gf2_accel16_sse_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			    gf2_matrix *dest);

/* Calculate src*mm from a gf2_mulmat acceleration structure for mm.
** Destination CAN NOT be the same as source.
*/
gf2_matrix *gf2_matrix_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			      gf2_matrix *dest);

void gf2_accel4_x86_mulmat(const gf2_matrix *src,const gf2_mulmat *mm,
			   gf2_matrix *dest);


/* Free a gf2_mulmat structure */
void gf2_free_mulmat(gf2_mulmat *mm);


/* Calculate an orthonormal representation for the vector space represented
** by the matrix a. Put the result in d, or allocate a new matrix for it.
** Return d or the new matrix. 
**
** Orthonormal representation has unit vectors 
*/
gf2_matrix *gf2_matrix_orthonormalize(const gf2_matrix *a,
				      gf2_matrix *d);

#endif
