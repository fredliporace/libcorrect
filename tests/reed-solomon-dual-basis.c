#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

#include "../include/correct/reed-solomon/reed-solomon.h"
#include "rs_tester.h"

void print_test_type(size_t block_length, size_t message_length,
                     size_t num_errors, size_t num_erasures) {
    printf(
        "testing reed solomon block length=%zu, message length=%zu, "
        "errors=%zu, erasures=%zu...",
        block_length, message_length, num_errors, num_erasures);
}

void fail_test() {
    printf("FAILED\n");
    exit(1);
}

void pass_test() { printf("PASSED\n"); }

void run_tests(correct_reed_solomon *rs, rs_testbench *testbench,
               size_t block_length, size_t test_msg_length, size_t num_errors,
               size_t num_erasures, size_t num_iterations) {
    rs_test test;
    test.encode = rs_correct_encode;
    test.decode = rs_correct_decode;
    test.encoder = rs;
    test.decoder = rs;
    print_test_type(block_length, test_msg_length, num_errors, num_erasures);
    for (size_t i = 0; i < num_iterations; i++) {
        rs_test_run run = test_rs_errors(&test, testbench, test_msg_length,
                                         num_errors,
                                     num_erasures);
        if (!run.output_matches) {
            fail_test();
        }
    }
    pass_test();
}

/*
  Generates dual basis representation lookuptables.
  Adapted from https://raw.githubusercontent.com/Opendigitalradio/
  ka9q-fec/master/gen_ccsds_tal.c
 */
/* Alpha to dual */
uint8_t Taltab[256];
/* Dual to alpha */
uint8_t Tal1tab[256];

void generate_dual_basis_lookuptables(void) {

  static uint8_t tal[] = { 0x8d, 0xef, 0xec, 0x86, 0xfa, 0x99, 0xaf, 0x7b };

  /* Generate conversion lookup tables between conventional alpha representation
   * (@**7, @**6, ...@**0)
   *  and Berlekamp's dual basis representation
   * (l0, l1, ...l7)
   */
  int i,j,k;

  for(i=0;i<256;i++){/* For each value of input */
    Taltab[i] = 0;
    for(j=0;j<8;j++) /* for each column of matrix */
      for(k=0;k<8;k++){ /* for each row of matrix */
	if(i & (1<<k))
	   Taltab[i] ^= tal[7-k] & (1<<j);
      }
    Tal1tab[Taltab[i]] = i;
  }
  printf("unsigned char Taltab[] = {\n");
  for(i=0;i<256;i++){
    if((i % 16) == 0)
      printf("\n");
    printf("0x%02x,",Taltab[i]);
  }
  printf("\n};\n\nunsigned char Tal1tab[] = {");
  for(i=0;i<256;i++){
    if((i % 16) == 0)
      printf("\n");
    printf("0x%02x,",Tal1tab[i]);
  }
  printf("\n};\n");
}

void show_block(uint8_t* block) {
  for(size_t i = 0 ; i < 255; i++) {
    printf(" %02X", block[i]);
  }
}

void dual_to_alpha(uint8_t* dual, uint8_t* alpha) {
  for(size_t i = 0; i < 255; i++) {
    alpha[i] = Tal1tab[dual[i]];
  }
}

/*
  Reads a CCSDS 1028 byte packet from file and split its content in
  4 blocks (interleave 4).
 */
uint8_t** read_blocks_from_file( char* input_filename ) {
  /* Ptrs to blocks */
  uint8_t** block_ptrs = calloc(4, sizeof(uint8_t*));
  for(size_t block = 0; block < 4; block++) {
    block_ptrs[block] = calloc(255, 1);
  }

  /* Read CCSDS packet from file. 4 bytes sync + 1024 bytes */
  FILE* fp = fopen(input_filename, "rb");
  assert( fp );
  uint8_t* ccsds_packet = calloc(1028, 1);
  size_t bytes_read = fread(ccsds_packet, 1, 1028, fp);
  assert( bytes_read == 1028 );

  size_t index = 4; // Skip the sync word
  size_t block_index = 0;
  while(index < 1028) {
    for(size_t block = 0; block < 4; block++) {
      block_ptrs[block][block_index] = ccsds_packet[index++];
    }
    block_index++;
  }

  fclose(fp);
  free(ccsds_packet);
  return block_ptrs;
}

int main() {
    srand(time(NULL));

    generate_dual_basis_lookuptables();
    uint8_t** blocks = read_blocks_from_file("ccsds_packets.dat");

    size_t block_length = 255;
    size_t min_distance = 16;
    size_t message_length = block_length - min_distance;

    correct_reed_solomon *rs = correct_reed_solomon_create(
        correct_rs_primitive_polynomial_ccsds, 120, 11, min_distance);

    printf("RS Generator poly coefficients\n");
    for(size_t i = 0; i < rs->min_distance + 1; i++ ) {
      printf("%02d ", rs->generator.coeff[i]);
    }
    printf("\n");

    printf("Original block\n");
    show_block(blocks[0]);
    printf("\n");

    uint8_t work_block[256];
    uint8_t reference_block[256];
    dual_to_alpha(blocks[0], work_block);
    dual_to_alpha(blocks[0], reference_block);

    // Include errors in work block
    work_block[3] = 0x5e;
    work_block[7] = 0x5e;
    work_block[9] = 0x5e;
    work_block[16] = 0x5e;
    work_block[200] = 0x5e;
    work_block[101] = 0x5e;
    work_block[43] = 0x5e;
    work_block[34] = 0x5e;
    //work_block[35] = 0x5e;

    printf("Dual block\n");
    show_block(work_block);
    printf("\n");

    uint8_t checksymbols[256];
    ssize_t block_decoded;
    block_decoded = correct_reed_solomon_decode(rs, work_block, block_length,
                                                checksymbols);
    printf("block_decoded: %d\n", block_decoded);
    assert(block_decoded == message_length);

    printf("Checksymbols (original, computed)\n");
    for(size_t i = 0; i < message_length; i++) {
      printf("\t%02X %02X\n", work_block[i], checksymbols[i]);
      assert(reference_block[i] == checksymbols[i]);
    }
    printf("\n");

    printf("test passed\n");
    return 0;

    exit(0);

    rs_testbench *testbench = rs_testbench_create(block_length, min_distance);

    run_tests(rs, testbench, block_length, message_length / 2, 0, 0, 20000);
    run_tests(rs, testbench, block_length, message_length, 0, 0, 20000);
    run_tests(rs, testbench, block_length, message_length / 2, min_distance / 2,
              0, 20000);
    run_tests(rs, testbench, block_length, message_length, min_distance / 2, 0,
              20000);
    run_tests(rs, testbench, block_length, message_length / 2, 0, min_distance,
              20000);
    run_tests(rs, testbench, block_length, message_length, 0, min_distance,
              20000);
    run_tests(rs, testbench, block_length, message_length / 2, min_distance / 4,
              min_distance / 2, 20000);
    run_tests(rs, testbench, block_length, message_length, min_distance / 4,
              min_distance / 2, 20000);

    rs_testbench_destroy(testbench);
    correct_reed_solomon_destroy(rs);

}
