#include <flame.h>
#include <stdio.h>
#include <random>

void printArray( int* arr, int n )
{
    printf( "[ " );
    for( int i = 0; i < n; i++)
        printf( "%i, ", arr[i] );
    printf( "]" );
}

int main()
{
  float **data = nullptr;
  int i, j, N = 1000, M = 2;
  FILE *fin;

//  fin = fopen( "C:\\Users\\Lazar\\Desktop\\flame-clustering-master\\matrix.txt", "r" );
//  fscanf( fin, "%i %i\n", & N, & M );

  //  printf( "Reading dataset with %i rows and %i columns\n", N, M );

  printf( "Creating dataset with %i rows and %i columns\n", N, M );

  data = new float*[ N ];
  for( i=0; i<N; i++ ){
    data[i] = new float[ M ];
    for( j=0; j<M; j++ )
        data[i][j] = (rand() % 100000) / 100.0;
//        fscanf( fin, "%f ", & data[i][j] );
//    printf("data: %f, %f \n", data[i][0], data[i][1] );
  }

  printf( "Putting data into Flame class.\n" );
  Flame flame = Flame( data, N, M );

  printf( "Detecting Cluster Supporting Objects ..." );
  fflush(stdout);
  flame.defineSupports( 10, .5 );
  printf( "done, found %i\n", flame.getCsoCount() );

  printf( "Propagating fuzzy memberships ... " );
  fflush( stdout );
  flame.localApproximation( 500, 1e-6 );
  printf( "done\n" );


  printf( "Defining clusters from fuzzy memberships ... " );
  fflush( stdout );
  flame.makeClusters( -1.0 );
  printf( "done\n" );

  printArray( flame.getClusters(), flame.getN() );



  printf( "All done.\n" );

  return 0;
}
