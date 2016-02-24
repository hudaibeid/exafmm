#ifndef HILBERT_C
#define HILBERT_C
//+++++++++++++++++++++++++++ PUBLIC-DOMAIN SOFTWARE ++++++++++++++++++++++++++
// Functions: TransposetoAxes AxestoTranspose
// Purpose: Transform in-place between Hilbert transpose and geometrical axes
// Example: b=5 bits for each of n=3 coordinates.
// 15-bit Hilbert integer = A B C D E F G H I J K L M N O is stored
// as its Transpose
// X[0] = A D G J M X[2]|
// X[1] = B E H K N <-------> | /X[1]
// X[2] = C F I L O axes |/
// high low 0------ X[0]
// Axes are stored conventially as b-bit integers.
// Author: John Skilling 20 Apr 2001 to 11 Oct 2003
//-----------------------------------------------------------------------------
//typedef uint32_t coord_t; // char,short,int for up to 8,16,32 bits per word


template<typename ArrayType>
void TransposetoAxes(ArrayType& X, int b, int n) // position, #bits, dimension
{ 
	typedef typename ArrayType::value_type coord_t;
	coord_t N = 2 << (b-1), P, Q, t;
	int i;
	// Gray decode by H ^ (H/2)
	t = X[n-1] >> 1;
	for( i = n-1; i > 0; i-- ) 
		X[i] ^= X[i-1];
	X[0] ^= t;
	// Undo excess work
	for( Q = 2; Q != N; Q <<= 1 ) {
		P = Q - 1;
		for( i = n-1; i >= 0 ; i-- )
			if( X[i] & Q ) 
				X[0] ^= P; // invert
			else
			{ 
				t = (X[0]^X[i]) & P;
				X[0] ^= t; 
				X[i] ^= t; 
			} // exchange
	} 
}

template<typename ArrayType>
void AxestoTranspose(ArrayType& X, int order, int dim)  // position, #bits, dimension
{ 
	typedef typename ArrayType::value_type coord_t;
	coord_t M = 1 << (order-1), P, Q, t;	
	// Inverse undo
	for( Q = M; Q > 1; Q >>= 1 ) {
		P = Q - 1;
		for(int i = 0; i < dim; i++ )
			if( X[i] & Q ) 
				X[0] ^= P; 																				// invert			
			else { 																					  	// exchange
				t = (X[0]^X[i]) & P; 
				X[0] ^= t; 
				X[i] ^= t; 
			} 
		} 
		
		// Gray encode
		for(int i = 1; i < dim; i++ ) 
			X[i] ^= X[i-1];
		
		t = 0;
		
		for( Q = M; Q > 1; Q >>= 1 )
			if( X[dim-1] & Q ) 
				t ^= Q-1;

		for(int i = 0; i < dim; i++ ) 
			X[i] ^= t;			
}




#endif
