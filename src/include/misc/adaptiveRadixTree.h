#ifndef MISC_ADAPTIVE_RADIX_TREE_H
#define MISC_ADAPTIVE_RADIX_TREE_H

// Copyright (c) 2012, Armon Dadgar
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of the organization nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL ARMON DADGAR BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Portability note: this was originally designed around the assumption that chars
// are 8 bits, and hence I've forced it to take uint8_ts instead. If this were to
// change, one could go through in the innards and use different node sizes,
// particularly with the node256, which uses the uint8_t values to directly index
// into an array.


#include "stddef.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct misc_art_node misc_art_node;

typedef struct {
  misc_art_node* root;
  misc_size_t size;
} misc_art_tree;

void misc_art_initialize(misc_art_tree* t);
int misc_art_invalidate(misc_art_tree* t);
misc_art_tree* misc_art_create();
int misc_art_destroy(misc_art_tree* t);


#if defined(__GNUC__) && !defined(__clang__) && __STDC_VERSION__ >= 199901L && 402 == (__GNUC__ * 100 + __GNUC_MINOR__)
#  define BROKEN_GCC_C99_INLINE
#  define misc_art_getSize(t) ((t)->size)
#else
inline misc_size_t misc_art_getSize(const misc_art_tree* t) {
  return t->size;
}
#endif

// returns the value that got overwritten, if any exists; otherwise NULL
// if NULL, errno will have been set, check for error
void* misc_art_insert(misc_art_tree* restrict t, const uint8_t* restrict key, misc_size_t keyLength, const void* restrict value);
// returns the value that got deleted, if any exists; otherwise NULL
// if NULL, errno will have been set, check for error
void* misc_art_delete(misc_art_tree* restrict t, const uint8_t* restrict key, misc_size_t keyLength);
// returns the value found, if any exists; otherwise NULL
// if NULL, errno will have been set, check for error
void* misc_art_search(const misc_art_tree* restrict t, const uint8_t* restrict key, misc_size_t keyLength);

typedef int (*misc_art_callback)(void* restrict data, const uint8_t* restrict key, misc_size_t keyLength, void* restrict value);
int misc_art_map(const misc_art_tree* restrict t, misc_art_callback cb, void* restrict data);
int misc_art_mapOverPrefix(const misc_art_tree* restrict t, const uint8_t* prefix, misc_size_t prefixLength, misc_art_callback cb, void* restrict data);

void misc_art_print(const misc_art_tree* t);

#ifdef __cplusplus
}
#endif

#endif // MISC_ADAPTIVE_RADIX_TREE_H

