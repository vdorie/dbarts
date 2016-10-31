#ifndef EXT_ADAPTIVE_RADIX_TREE_H
#define EXT_ADAPTIVE_RADIX_TREE_H

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

typedef struct ext_art_node ext_art_node;

typedef struct {
  ext_art_node* root;
  ext_size_t size;
} ext_art_tree;

void ext_art_initialize(ext_art_tree* t);
int ext_art_invalidate(ext_art_tree* t);
ext_art_tree* ext_art_create();
int ext_art_destroy(ext_art_tree* t);


#if defined(__GNUC__) && !defined(__clang__) && __STDC_VERSION__ >= 199901L && 402 == (__GNUC__ * 100 + __GNUC_MINOR__)
#  define BROKEN_GCC_C99_INLINE
#  define ext_art_getSize(t) ((t)->size)
#else
inline ext_size_t ext_art_getSize(const ext_art_tree* t) {
  return t->size;
}
#endif

// returns the value that got overwritten, if any exists; otherwise NULL
// if NULL, errno may have been set
void* ext_art_insert(ext_art_tree* restrict t, const uint8_t* restrict key, ext_size_t keyLength, const void* restrict value);
// returns the value that got deleted, if any exists; otherwise NULL
// if NULL, errno may have been set
void* ext_art_delete(ext_art_tree* restrict t, const uint8_t* restrict key, ext_size_t keyLength);
// returns the value found, if any exists; otherwise NULL
// if NULL, errno may have been set
void* ext_art_search(const ext_art_tree* restrict t, const uint8_t* restrict key, ext_size_t keyLength);

typedef int (*ext_art_callback)(void* restrict data, const uint8_t* restrict key, ext_size_t keyLength, void* restrict value);
int ext_art_map(const ext_art_tree* restrict t, ext_art_callback cb, void* restrict data);
int ext_art_mapOverPrefix(const ext_art_tree* restrict t, const uint8_t* prefix, ext_size_t prefixLength, ext_art_callback cb, void* restrict data);

void ext_art_print(const ext_art_tree* t);

#ifdef __cplusplus
}
#endif

#endif // EXT_ADAPTIVE_RADIX_TREE_H

