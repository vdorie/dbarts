#include <external/adaptiveRadixTree.h>
#include "config.h"

#include <errno.h>
#include <stdbool.h>
#include <external/stddef.h>
#include <stdlib.h> // malloc, free
#include <string.h> // memcpy, memcmp
#ifdef HAVE_SSE2
#  ifdef _MSC_VER
#    include <intrin.h>
#  else
#    include <emmintrin.h> // SSE2 intrinsics
#  endif
#endif

#include <external/io.h>

#if defined(__GNUC__) && (\
  (!defined(__clang__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6))) || \
  ( defined(__clang__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 7))))
#  define SUPPRESS_DIAGNOSTIC 1
#endif

// chosen to cause the keys to align to a 16 byte address, assuming
// that the base does as well
#if SIZEOF_SIZE_T == 8
#  define MAX_PARTIAL_LENGTH ((size_t) 22)
# else 
#  define MAX_PARTIAL_LENGTH ((size_t) 10)
#endif

// we define partial as part of the prefix, prefix being
// the whole shared initial string; so partial is the prefix
// if the prefix has length <= MAX_PARTIAL_LENGTH
struct ext_art_node {
  uint8_t type;
  uint8_t numChildren;
  uint8_t partial[MAX_PARTIAL_LENGTH];
  size_t prefixLength;
};

typedef ext_art_tree Tree;
typedef ext_art_node Node;

typedef enum {
  NODE4 = 0,
  NODE16,
  NODE48,
  NODE256
} NodeType;

typedef struct {
  Node n;
  uint8_t keys[4];
  Node* children[4];
} Node4;

typedef struct {
  Node n;
  uint8_t keys[16];
  Node* children[16];
} Node16;

typedef struct {
  Node n;
  uint8_t keys[256];
  Node* children[48];
} Node48;

typedef struct {
  Node n;
  Node* children[256];
} Node256;

typedef struct {
  const void* value;
  size_t keyLength;
  const uint8_t key[];
} Leaf;

#if ALIGNOF_VOIDP > 0
#  define TAG_BIT 1
#elif SIZEOF_SIZE_T
// this is in no way guaranteed to work, and maybe there is a way to query
// the base and direction of allocation
#  define TAG_BIT 0x8000000000000000ull
#else
#  error Cannot find a safe place to put pointer tag
#endif

#define nodeIsLeaf(_X_) ((uintptr_t) _X_ & TAG_BIT)
#define tagNodeAsLeaf(_X_)     ((void*) ((uintptr_t) _X_ |  TAG_BIT))
#define getRawLeafPointer(_X_) ((void*) ((uintptr_t) _X_ & ~TAG_BIT))

#define INVALID_LENGTH ((size_t) -1)

// forward declarations
static int insertKeyValue(Node* restrict n, const uint8_t* restrict key, size_t keyLength, const void* restrict value, size_t depth,
                          Node* restrict* restrict positionInParent, bool* restrict addedNode, void* restrict* restrict oldValue);
static int deleteKey(Node* restrict n, const uint8_t* restrict key, size_t keyLength, size_t depth,
                     Node* restrict* restrict positionInParent, Leaf* restrict* restrict deletedLeaf);

static bool keysMatch(const Leaf* restrict l, const uint8_t* restrict key, size_t keyLength);
static bool prefixesMatch(const Leaf* restrict l, const uint8_t* restrict prefix, size_t prefixLength);
static size_t getLongestCommonPrefixLength(const Leaf* restrict l1, const Leaf* restrict l2, size_t depth);
// first one only uses node's partial, second one digs down tree to find full prefix
static size_t getPartialMismatchIndex(const Node* restrict n, const uint8_t* restrict key, size_t keyLength, size_t depth);
static size_t getPrefixMismatchIndex(const Node* restrict n, const uint8_t* restrict key, size_t keyLength, size_t depth);

static Node** findChildMatchingKey(const Node* n, uint8_t c);
static Leaf* getMinimumLeafUnderNode(const Node* n);
// static Leaf* getMaximumLeafUnderNode(const Node* n);

static int map(const Node* restrict n, ext_art_callback cb, void* restrict data);

static int addChild(Node* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent);
static int removeChild(Node* restrict n, uint8_t c, Node* restrict* restrict positionInParent, Node* restrict* restrict leafRef);
static int addChild4(Node4* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent);

static int destroyNode(Node* n);
static Node* createNode(NodeType type);
static Leaf* createLeaf(const uint8_t* restrict key, size_t keyLength, const void* restrict value);
static void copyNodeHeader(Node* restrict dest, const Node* restrict src);

static void printAtDepth(Node* const n, size_t depth);

static inline size_t getMinLength(size_t a, size_t b) { return (a < b) ? a : b; }

#ifdef __GCC__
#  define countTrailingZeros(_X_) __builtin_ctz(_X_)
#elif defined(_MSC_VER)
// it would be possible to use __lzcnt when ABM is available, but that is relatively recent
#  include <intrin.h>
static unsigned long __inline countTrailingZeros(unsigned long x) {
  unsigned long result;
  _BitScanForward(&result, x);
  return result;
}
#elif defined(__INTEL_COMPILER)
#  include <immintrin.h>
#  define countTrailingZeros(_X_) _bit_scan_forward(_X_)
#elif defined(HAVE_FFS)
#  include <strings.h>
#  define countTrailingZeros(_X_) (ffs(_X_) - 1)
#else
static inline uint32_t countTrailingZeros(uint32_t x) {
  static const uint8_t ctzTable[] = {
    0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 
    0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 
    0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 
    0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 
    0, 1, 0, 2, 0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 
    0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 
    0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 
    0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1 };

  uint32_t n = 0;
  if ((x & 0x0000FFFF) == 0) { n += 16; x >>= 16; }
  if ((x & 0x000000FF) == 0) { n +=  8; x >>= 8;  }
  n += ctzTable[x];
  return n;
#endif

void ext_art_initialize(Tree* t) {
  t->root = NULL;
  t->size = 0;
}

int ext_art_invalidate(Tree* t) {
  int result = destroyNode(t->root);
  t->root = NULL;
  t->size = 0;
  return result;
}

Tree* ext_art_create() {
  Tree* result = (Tree*) malloc(sizeof(Tree));
  if (result != NULL) ext_art_initialize(result);
  return result;
}

int ext_art_destroy(Tree* t) {
  if (t == NULL) return 0;
  int result = ext_art_invalidate(t);
  free(t);
  return result;
}

void* ext_art_insert(Tree* restrict t, const uint8_t* restrict key, size_t keyLength, const void* restrict value)
{
  bool addedNode = true;
  void* oldValue = NULL;
  int errorCode = insertKeyValue(t->root, key, keyLength, value, 0, &t->root, &addedNode, &oldValue);
  if (errorCode != 0) { errno = errorCode; return NULL; }
  if (addedNode == true) t->size++;
  return oldValue;
}

void* ext_art_delete(ext_art_tree* restrict t, const uint8_t* restrict key, size_t keyLength)
{
  Leaf* l = NULL;
  int errorCode = deleteKey(t->root, key, keyLength, 0, &t->root, &l);
  if (errorCode != 0) { errno = errorCode; return NULL; }
  
  if (l == NULL) return NULL;
  
  t->size--;
  void* oldValue = (void*) l->value;
  free(l);
  return oldValue;
}

void* ext_art_search(const Tree* restrict t, const uint8_t* restrict key, size_t keyLength)
{
  Node** child;
  const Node* n = t->root;
  size_t depth = 0;
  while (n != NULL) {
    if (nodeIsLeaf(n)) {
      n = (const Node*) getRawLeafPointer(n);
      if (keysMatch((const Leaf*) n, key, keyLength))
        return (void*) ((Leaf*) n)->value;
      return NULL;
    }
    
    if (n->prefixLength > 0) {
      size_t partialLength = getPartialMismatchIndex(n, key, keyLength, depth);
      if (partialLength != getMinLength(MAX_PARTIAL_LENGTH, n->prefixLength))
        return NULL;
      depth = depth + n->prefixLength;
    }
    
    child = findChildMatchingKey(n, key[depth]);
    if (child != NULL) {
      n = (const Node*) *child;
    } else {
      n = NULL;
      // errno will propagate
    }
    ++depth;
  }
  return NULL;
}

int ext_art_map(const Tree* restrict t, ext_art_callback cb, void* restrict data) {
  return map(t->root, cb, data);
}

int ext_art_mapOverPrefix(const ext_art_tree* restrict t, const uint8_t* prefix, ext_size_t prefixLength, ext_art_callback cb, void* restrict data)
{
  size_t depth = 0;
  const Node* n = t->root;
  while (n != NULL) {
    if (nodeIsLeaf(n)) {
      const Leaf* l = (const Leaf*) getRawLeafPointer(n);
      if (prefixesMatch(l, prefix, prefixLength))
        return cb(data, l->key, l->keyLength, (void*) l->value);
      else
        return 0;
    }
    
    if (depth == prefixLength) {
      const Leaf* l = getMinimumLeafUnderNode(n);
      if (prefixesMatch(l, prefix, prefixLength))
        return map(n, cb, data);
      else
        return 0;
    }
    
    if (n->prefixLength > 0) {
      size_t commonPrefixLength = getPrefixMismatchIndex(n, prefix, prefixLength, depth);
      
      if (commonPrefixLength == 0 || commonPrefixLength == INVALID_LENGTH)
        return 0;
      else if (depth + commonPrefixLength == prefixLength)
        return map(n, cb, data);
      
      depth += n->prefixLength;
    }
    
    Node** child = findChildMatchingKey(n, prefix[depth]);
    n = (child != NULL) ? *child : NULL;
    ++depth;
  }
  return 0;
}

void ext_art_print(const Tree* t)
{
  printAtDepth(t->root, 0);
}



static int insertKeyValue(Node* restrict n, const uint8_t* restrict key, size_t keyLength, const void* restrict value, size_t depth,
                          Node* restrict* restrict positionInParent, bool* restrict addedNode, void* restrict* restrict oldValue)
{
  if (n == NULL) {
    Leaf* l = createLeaf(key, keyLength, value);
    if (l == NULL) return errno;
    *positionInParent = (Node*) tagNodeAsLeaf(l);
    return 0;
  }

  if (nodeIsLeaf(n)) {
    Leaf* l = (Leaf*) getRawLeafPointer(n);
    
    if (keysMatch(l, key, keyLength)) {
      *addedNode = false;
      *oldValue = (void*) l->value;
      l->value = value;
      return 0;
    }
    
    // recursed to a leaf node that differs from given key, replace leaf by internal node
    // and point to two new leaves
    Node4* newNode = (Node4*) createNode(NODE4);
    if (newNode == NULL) return errno;
    
    Leaf* l2 = createLeaf(key, keyLength, value);
    if (l2 == NULL) return destroyNode((Node*) newNode), errno;
    
    size_t prefixLength = getLongestCommonPrefixLength(l, l2, depth);
    newNode->n.prefixLength = prefixLength;
    memcpy(newNode->n.partial, key + depth, getMinLength(prefixLength, MAX_PARTIAL_LENGTH));
    
    int errorCode = 0;
    if ((errorCode = addChild4(newNode, l->key[depth + prefixLength],  tagNodeAsLeaf(l), NULL)) != 0 ||
        (errorCode = addChild4(newNode, l2->key[depth + prefixLength], tagNodeAsLeaf(l2), NULL) != 0))
    {
      free(l2);
      destroyNode((Node*) newNode);
      return errorCode;
    }
    *positionInParent = (Node*) newNode;
    return 0;
  }
  
  // check if node has a prefix
  if (n->prefixLength > 0) {
    // see if prefixes differ up to this point: if so, split; if not, recurse
    size_t commonPrefixLength = getPrefixMismatchIndex(n, key, keyLength, depth);
    if (commonPrefixLength == INVALID_LENGTH) return EINVAL;
    
    if (commonPrefixLength < n->prefixLength) {
      Node4* newNode = (Node4*) createNode(NODE4);
      if (newNode == NULL) return errno;
      newNode->n.prefixLength = commonPrefixLength;
      memcpy(newNode->n.partial, n->partial, getMinLength(MAX_PARTIAL_LENGTH, commonPrefixLength));
      
      // adjust the prefix of the old node
      if (n->prefixLength <= MAX_PARTIAL_LENGTH) {
        if (addChild4(newNode, n->partial[commonPrefixLength], n, NULL) != 0) return destroyNode((Node*) newNode), errno;
        n->prefixLength -= (commonPrefixLength + 1);
        memmove(n->partial, n->partial + commonPrefixLength + 1, getMinLength(MAX_PARTIAL_LENGTH, n->prefixLength));
      } else {
        n->prefixLength -= (commonPrefixLength + 1);
        Leaf* l = getMinimumLeafUnderNode(n);
        if (addChild4(newNode, l->key[depth + commonPrefixLength], n, NULL) != 0) return destroyNode((Node*) newNode), errno;
        memcpy(n->partial, l->key + depth + commonPrefixLength + 1, getMinLength(MAX_PARTIAL_LENGTH, n->prefixLength));
      }
      
      // insert the new leaf
      Leaf* l = createLeaf(key, keyLength, value);
      if (l == NULL) return errno;
      int errorCode = 0;
      if ((errorCode = addChild4(newNode, key[depth + commonPrefixLength], tagNodeAsLeaf(l), NULL)) != 0) {
        free(l);
        return destroyNode((Node*) newNode), errorCode;
      }
      
      *positionInParent = (Node*) newNode;
      return 0;
    } else {
     depth += n->prefixLength;
    }
  }

  errno = 0;
  Node** child = findChildMatchingKey(n, key[depth]);
  if (child != NULL)
    return insertKeyValue(*child, key, keyLength, value, depth + 1, child, addedNode, oldValue);
  
  if (errno != 0) return errno;

  // no child found, node goes within us
  Leaf* l = createLeaf(key, keyLength, value);
  if (l == NULL) return errno;
  return addChild(n, key[depth], tagNodeAsLeaf(l), positionInParent);
}

static int deleteKey(Node* restrict n, const uint8_t* restrict key, size_t keyLength, size_t depth,
                     Node* restrict* restrict positionInParent, Leaf* restrict* restrict deletedLeaf)
{
  if (n == NULL) return 0;
  
  if (nodeIsLeaf(n)) {
    Leaf* l = (Leaf*) getRawLeafPointer(n);
    if (keysMatch(l, key, keyLength)) {
      *positionInParent = NULL;
      *deletedLeaf = l;
      return 0;
    }
  }
  
  // scan the partial to see if we can abort quickly, without checking leaves
  if (n->prefixLength > 0) {
    size_t commonPartialLength = getPartialMismatchIndex(n, key, keyLength, depth);
    if (commonPartialLength != getMinLength(MAX_PARTIAL_LENGTH, n->prefixLength)) return 0;
    depth += n->prefixLength;
  }
  
  errno = 0;
  Node** child = findChildMatchingKey(n, key[depth]);
  if (child == NULL) return errno;
  
  if (nodeIsLeaf(*child)) {
    Leaf* l = (Leaf*) getRawLeafPointer(*child);
    if (keysMatch(l, key, keyLength)) {
      int errorCode = removeChild(n, key[depth], positionInParent, child);
      if (errorCode != 0) return errorCode;
      
      *deletedLeaf = l;
      return 0;
    }
    return 0;
  }
  
  return deleteKey(*child, key, keyLength, depth + 1, child, deletedLeaf);
}

static bool keysMatch(const Leaf* restrict l, const uint8_t* restrict key, size_t keyLength) {
  if (l->keyLength != keyLength) return false;

  return memcmp(l->key, key, keyLength) == 0;
}

static bool prefixesMatch(const Leaf* restrict l, const uint8_t* restrict prefix, size_t prefixLength) {
  if (l->keyLength < prefixLength) return false;
  
  return memcmp(l->key, prefix, prefixLength) == 0;
}

static size_t getPartialMismatchIndex(const Node* restrict n, const uint8_t* restrict key, size_t keyLength, size_t depth) {
  size_t maxIndex = getMinLength(getMinLength(MAX_PARTIAL_LENGTH, n->prefixLength), keyLength - depth);
  size_t index;
  for (index = 0; index < maxIndex; ++index) {
    if (n->partial[index] != key[depth + index]) return index;
  }
  return index;
}

static size_t getLongestCommonPrefixLength(const Leaf* restrict l1, const Leaf* restrict l2, size_t depth) {
  size_t maxIndex = getMinLength(l1->keyLength, l2->keyLength) - depth;
  size_t index;
  for (index = 0; index < maxIndex; ++index) {
    if (l1->key[depth + index] != l2->key[depth + index])
      return index;
  }
  return index;
}

static size_t getPrefixMismatchIndex(const Node* restrict n, const uint8_t* restrict key, size_t keyLength, size_t depth) {
  size_t maxIndex = getMinLength(getMinLength(MAX_PARTIAL_LENGTH, n->prefixLength), keyLength - depth);
  size_t index;
  for (index = 0; index < maxIndex; ++index) {
    if (n->partial[index] != key[depth + index])
      return index;
  }

  // If the prefix is short we can avoid finding a leaf
  if (n->prefixLength <= MAX_PARTIAL_LENGTH) return index;
  
  // Prefix is longer than what we've checked, find a leaf
  Leaf* l = getMinimumLeafUnderNode(n);
  if (l == NULL) return INVALID_LENGTH;
  
  maxIndex = getMinLength(l->keyLength, keyLength) - depth;
  for (; index < maxIndex; ++index) {
    if (l->key[index + depth] != key[depth + index])
      return index;
  }
  
  return index;
}

static Node** findChildMatchingKey(const Node* n, uint8_t c)
{
  uint8_t index;
  union {
    const Node4  * p1;
    const Node16 * p2;
    const Node48 * p3;
    const Node256* p4;
  } p;
  switch (n->type) {
    case NODE4:
    p.p1 = (const Node4*) n;
    for (index = 0; index < n->numChildren; ++index)
      if (p.p1->keys[index] == c) return (Node**) &p.p1->children[index];
    break;
    
    case NODE16:
    p.p2 = (const Node16*) n;
    if (((uintptr_t) p.p2->keys % 0x10) != 0) ext_throwError("BAAAALS!!!!\n");
#ifdef HAVE_SSE2
    {
      __m128i comparison = _mm_cmpeq_epi8(_mm_set1_epi8(c), _mm_loadu_si128((__m128i*) p.p2->keys));
      
      unsigned int bitfield = ((1 << n->numChildren) - 1) & (unsigned int) _mm_movemask_epi8(comparison);
      if (bitfield != 0) return (Node**) &p.p2->children[countTrailingZeros(bitfield)];
    }
#else
    p.p2 = (const Node16*) n;
    for (index = 0; index < n->numChildren; ++index)
      if (p.p2->keys[index] == c) return (Node**) &p.p2->children[index];
#endif
    break;
    
    case NODE48:
    p.p3 = (const Node48*) n;
    index = p.p3->keys[c];
    if (index != 0) return (Node**) &p.p3->children[index - 1];
    break;
    
    case NODE256:
    p.p4 = (const Node256*) n;
    if (p.p4->children[c] != NULL) return (Node**) &p.p4->children[c];
    break;
    
    default:
    errno = EINVAL;
    break;
  }
  return NULL;
}

static Leaf* getMinimumLeafUnderNode(const Node* n) {
  if (n == NULL) return NULL;
  if (nodeIsLeaf(n)) return getRawLeafPointer(n);
  if (n->numChildren == 0) { errno = EINVAL; return NULL; }
  
  uint8_t index;
  switch (n->type) {
    case NODE4:
    return getMinimumLeafUnderNode(((Node4*) n)->children[0]);
    
    case NODE16:
    return getMinimumLeafUnderNode(((Node16*) n)->children[0]);
    
    case NODE48:
    index = 0;
    while (((Node48*) n)->keys[index] != 0) ++index;
    index = ((Node48*) n)->keys[index] - 1;
    return getMinimumLeafUnderNode(((Node48*) n)->children[index]);
    
    case NODE256:
    index = 0;
    while (((Node256*) n)->children[index] != 0) ++index;
    return getMinimumLeafUnderNode(((Node256*) n)->children[index]);
    
    default:
    errno = EINVAL;
    return NULL;
  }
}

/* static Leaf* getMaximumLeafUnderNode(const Node* n) {
  if (n == NULL) return NULL;
  if (nodeIsLeaf(n)) return getRawLeafPointer(n);
  if (n->numChildren == 0) { errno = EINVAL; return NULL; }
  
  uint8_t index;
  switch (n->type) {
    case NODE4:
    return getMaximumLeafUnderNode(((Node4*) n)->children[n->numChildren - 1]);
    
    case NODE16:
    return getMaximumLeafUnderNode(((Node16*) n)->children[n->numChildren - 1]);
    
    case NODE48:
    index = 255;
    while (((Node48*) n)->keys[index] != 0) --index;
    index = ((Node48*) n)->keys[index] - 1;
    return getMaximumLeafUnderNode(((Node48*) n)->children[index]);
    
    case NODE256:
    index = 255;
    while (((Node256*) n)->children[index] != 0) --index;
    return getMaximumLeafUnderNode(((Node256*) n)->children[index]);
    
    default:
    errno = EINVAL;
    return NULL;
  }
} */

static int map(const Node* restrict n, ext_art_callback cb, void* restrict data)
{
  if (n == NULL) return 0;
  if (nodeIsLeaf(n)) {
    const Leaf* l = (Leaf*) getRawLeafPointer(n);
    return cb(data, l->key, l->keyLength, (void*) l->value);
  }
  
  union {
    const Node4  * p1;
    const Node16 * p2;
    const Node48 * p3;
    const Node256* p4;
  } p;
  uint8_t index;
  int result;
  switch (n->type) {
    case NODE4:
    p.p1 = (const Node4*) n;
    for (index = 0; index < n->numChildren; ++index) {
      if (p.p1->keys[index] != 0 &&
          (result = map(p.p1->children[index], cb, data)) != 0)
        return result;
    }
    break;
    
    case NODE16:
    p.p2 = (const Node16*) n;
    for (index = 0; index < n->numChildren; ++index) {
      if (p.p2->keys[index] != 0 &&
          (result = map(p.p2->children[index], cb, data)) != 0)
        return result;
    }
    break;
    
    case NODE48:
    p.p3 = (const Node48*) n;
    for (size_t i = 0; i < 256; ++i) {
      index = p.p3->keys[i];
      if (index != 0) continue;
      
      if ((result = map(p.p3->children[index - 1], cb, data)) != 0)
        return result;
    }
    break;
    
    case NODE256:
    p.p4 = (const Node256*) n;
    for (size_t i = 0; i < 256; ++i) {
      if (p.p4->children[i] == NULL) continue;
      
      if ((result = map(p.p4->children[i], cb, data)) != 0)
        return result;
    }
    break;
  }
  return EINVAL;
}

static int addChild16(Node16* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent);
static int addChild48(Node48* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent);
static int addChild256(Node256* restrict n, uint8_t c, void* restrict child);

static int addChild(Node* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent)
{
  switch (n->type) {
    case NODE4:
    return addChild4((Node4*) n, c, child, positionInParent);
    case NODE16:
    return addChild16((Node16*) n, c, child, positionInParent);
    case NODE48:
    return addChild48((Node48*) n, c, child, positionInParent);
    case NODE256:
    return addChild256((Node256*) n, c, child);
  }
  return EINVAL;
}

static int addChild4(Node4* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent)
{
  if (n->n.numChildren < 4) {
    size_t index;
    for (index = 0; index < n->n.numChildren; ++index)
      if (c < n->keys[index]) break;
    
    // shift over old
    if (index < n->n.numChildren) {
      memmove(n->keys + index + 1,     n->keys + index,     (n->n.numChildren - index) * sizeof(uint8_t));
      memmove(n->children + index + 1, n->children + index, (n->n.numChildren - index) * sizeof(void*));
    }
    
    n->keys[index] = c;
    n->children[index] = child;
    n->n.numChildren++;
    
    return 0;
  }
  
  Node16* newNode = (Node16*) createNode(NODE16);
  if (newNode == NULL) return errno;
  
  memcpy(newNode->keys,     n->keys,     n->n.numChildren * sizeof(uint8_t));
  memcpy(newNode->children, n->children, n->n.numChildren * sizeof(void*));
  copyNodeHeader((Node*) newNode, (const Node*) n);
  *positionInParent = (Node*) newNode;
  free(n);
  
  return addChild16(newNode, c, child, NULL);
}

static int addChild16(Node16* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent)
{
  if (n->n.numChildren < 16) {
    size_t index;
#if defined(HAVE_SSE2)
    __m128i comparison = _mm_cmplt_epi8(_mm_set1_epi8(c), _mm_loadu_si128((__m128i*) n->keys));

    // Use a mask to ignore children that don't exist
    unsigned int mask = (1 << n->n.numChildren) - 1;
    unsigned int bitfield = _mm_movemask_epi8(comparison) & mask;

    // Check if less than any
    index = (bitfield != 0) ? countTrailingZeros(bitfield) : n->n.numChildren;
#else
    for (index = 0; index < n->n.numChildren; ++index)
      if (c < n->keys[index]) break;
#endif
    
    if (index < n->n.numChildren) {
      memmove(n->keys + index + 1,     n->keys + index,     (n->n.numChildren - index) * sizeof(uint8_t));
      memmove(n->children + index + 1, n->children + index, (n->n.numChildren - index) * sizeof(void*));
    }

    // Set the child
    n->keys[index] = c;
    n->children[index] = child;
    n->n.numChildren++;
    
    return 0;
  }
  
  Node48* newNode = (Node48*) createNode(NODE48);
  if (newNode == NULL) return errno;

  memcpy(newNode->children, n->children, n->n.numChildren * sizeof(void*));
  for (uint8_t i = 0; i < n->n.numChildren; ++i) {
    newNode->keys[n->keys[i]] = i + 1;
  }
  copyNodeHeader((Node*) newNode, (const Node*) n);
  *positionInParent = (Node*) newNode;
  free(n);
  return addChild48(newNode, c, child, NULL);
}

static int addChild48(Node48* restrict n, uint8_t c, void* restrict child, Node* restrict* restrict positionInParent)
{
  if (n->n.numChildren < 48) {
    size_t pos = 0;
    while (n->children[pos]) ++pos;
    n->children[pos] = child;
    n->keys[c] = pos + 1;
    n->n.numChildren++;
    
    return 0;
  }
  
  Node256* newNode = (Node256*) createNode(NODE256);
  if (newNode == NULL) return errno;
  
  for (size_t i = 0; i < 256; ++i) {
    if (n->keys[i] != 0) {
      newNode->children[i] = n->children[n->keys[i] - 1];
    }
  }
  copyNodeHeader((Node*) newNode, (const Node*) n);
  *positionInParent = (Node*) newNode;
  free(n);
  return addChild256(newNode, c, child);
}

static int addChild256(Node256* restrict n, uint8_t c, void* restrict child)
{
  n->n.numChildren++;
  n->children[c] = child;
  return 0;
}

static int removeChild4(Node4* restrict n, Node* restrict* restrict parentRef, Node* restrict* restrict leafRef);
static int removeChild16(Node16* restrict n, Node* restrict* restrict parentRef, Node* restrict* restrict leafRef);
static int removeChild48(Node48* restrict n, Node* restrict* restrict parentRef, uint8_t c);
static int removeChild256(Node256* restrict n, Node* restrict* restrict parentRef, uint8_t c);


static int removeChild(Node* restrict n, uint8_t c, Node* restrict* restrict parentRef, Node* restrict* restrict leafRef)
{
  switch (n->type) {
    case NODE4:
    return removeChild4((Node4*) n, parentRef, leafRef);
    case NODE16:
    return removeChild16((Node16*) n, parentRef, leafRef);
    case NODE48:
    return removeChild48((Node48*) n, parentRef, c);
    case NODE256:
    return removeChild256((Node256*) n, parentRef, c);
  }
  return EINVAL;
}

static int removeChild4(Node4* restrict n, Node* restrict* restrict parentRef, Node* restrict* restrict leafRef)
{
  size_t shift = leafRef - n->children;
  memmove(n->keys + shift,     n->keys + shift + 1    , (n->n.numChildren - 1 - shift) * sizeof(uint8_t));
  memmove(n->children + shift, n->children + shift + 1, (n->n.numChildren - 1 - shift) * sizeof(void*));
  n->n.numChildren--;
  
  if (n->n.numChildren == 1) {
    Node* child = n->children[0];
    if (!nodeIsLeaf(child)) {
      size_t prefixLength = n->n.prefixLength;
      if (prefixLength < MAX_PARTIAL_LENGTH) {
        n->n.partial[prefixLength] = n->keys[0];
        ++prefixLength;
      }
      if (prefixLength < MAX_PARTIAL_LENGTH) {
        size_t subPrefixLength = getMinLength(child->prefixLength, MAX_PARTIAL_LENGTH - prefixLength);
        memcpy(n->n.partial + prefixLength, child->partial, subPrefixLength);
        prefixLength += subPrefixLength;
      }
      
      memcpy(child->partial, n->n.partial, getMinLength(prefixLength, MAX_PARTIAL_LENGTH));
      child->prefixLength += n->n.prefixLength + 1;
    }
    *parentRef = child;
    free(n);
  }
  return 0;
}

static int removeChild16(Node16* restrict n, Node* restrict* restrict parentRef, Node* restrict* restrict leafRef)
{
  size_t shift = leafRef - n->children;
  memmove(n->keys + shift,     n->keys + shift + 1    , (n->n.numChildren - 1 - shift) * sizeof(uint8_t));
  memmove(n->children + shift, n->children + shift + 1, (n->n.numChildren - 1 - shift) * sizeof(void*));
  n->n.numChildren--;

  if (n->n.numChildren == 3) {
    Node4* newNode = (Node4*) createNode(NODE4);
    if (newNode == NULL) return errno;
    
    *parentRef = (Node*) newNode;
    copyNodeHeader((Node*) newNode, (const Node*) n);
    memcpy(newNode->keys,     n->keys,     4 * sizeof(uint8_t));
    memcpy(newNode->children, n->children, 4 * sizeof(void*));
    free(n);
  }
  return 0;
}

static int removeChild48(Node48* restrict n, Node* restrict* restrict parentRef, uint8_t c) {
  uint8_t pos = n->keys[c];
  n->keys[c] = 0;
  n->children[pos - 1] = NULL;
  n->n.numChildren--;

  if (n->n.numChildren == 12) {
    Node16* newNode = (Node16*) createNode(NODE16);
    if (newNode == NULL) return errno;
    
    *parentRef = (Node*) newNode;
    copyNodeHeader((Node*) newNode, (const Node*) n);

    uint8_t childIndex = 0;
    for (size_t i = 0; i < 256; ++i) {
      pos = n->keys[i];
      if (pos != 0) {
        newNode->keys[childIndex] = (uint8_t) i;
        newNode->children[childIndex] = n->children[pos - 1];
        ++childIndex;
      }
    }
    free(n);
  }
  return 0;
}

static int removeChild256(Node256* restrict n, Node* restrict* restrict parentRef, uint8_t c) {
  n->children[c] = NULL;
  n->n.numChildren--;

  if (n->n.numChildren == 37) {
    Node48* newNode = (Node48*) createNode(NODE48);
    if (newNode == NULL) return errno;
    
    *parentRef = (Node*) newNode;
    copyNodeHeader((Node*) newNode, (const Node*) n);

    uint8_t pos = 0;
    for (size_t i = 0; i < 256; ++i) {
      if (n->children[i] != NULL) {
        newNode->children[pos] = n->children[i];
        newNode->keys[i] = pos + 1;
        ++pos;
      }
    }
    free(n);
  }
  return 0;
}


static Node* createNode(NodeType type) {
  Node* n;
  switch (type) {
    case NODE4:
      n = (Node*) calloc(1, sizeof(Node4));
      break;
    case NODE16:
      n = (Node*) calloc(1, sizeof(Node16));
      break;
    case NODE48:
      n = (Node*) calloc(1, sizeof(Node48));
      break;
    case NODE256:
      n = (Node*) calloc(1, sizeof(Node256));
      break;
    default:
      return NULL;
  }
  if (n == NULL) return NULL;
  n->type = type;
  return n;
}

static int destroyNode(Node* n) {
  if (n == NULL) return 0;
  
  if (nodeIsLeaf(n)) {
    free(getRawLeafPointer(n));
    return 0;
  }
  
  union {
    Node4  * p1;
    Node16 * p2;
    Node48 * p3;
    Node256* p4;
  } p;
  switch (n->type) {
    case NODE4:
    p.p1 = (Node4*) n;
    for (uint8_t i = 0; i < n->numChildren; ++i) destroyNode(p.p1->children[i]);
    break;
    
    case NODE16:
    p.p2 = (Node16*) n;
    for (uint8_t i = 0; i < n->numChildren; ++i) destroyNode(p.p2->children[i]);
    break;
    
    case NODE48:
    p.p3 = (Node48*) n;
    for (uint8_t i = 0; i < n->numChildren; ++i) destroyNode(p.p3->children[i]);
    break;
    
    case NODE256:
    p.p4 = (Node256*) n;
    for (uint8_t i = 0; i < n->numChildren; ++i) destroyNode(p.p4->children[i]);
    break;
    
    default:
      return EINVAL;
  }
  
  free(n);
  
  return 0;
}

static Leaf* createLeaf(const uint8_t* restrict key, size_t keyLength, const void* restrict value) {
  Leaf* l = (Leaf*) malloc(sizeof(Leaf) + keyLength);
  if (l == NULL) return NULL;
  
  l->value = value;
  l->keyLength = keyLength;
  memcpy((uint8_t*) l->key, key, keyLength * sizeof(uint8_t));
  return l;
}

static void copyNodeHeader(Node* restrict dest, const Node* restrict src) {
  dest->numChildren = src->numChildren;
  dest->prefixLength = src->prefixLength;
  memcpy(dest->partial, src->partial, getMinLength(MAX_PARTIAL_LENGTH, src->prefixLength));
}

static void printAtDepth(Node* const n, size_t depth)
{
  if (n == NULL) {
    ext_printf("NULL\n");
    return;
  }
  
  if (nodeIsLeaf(n)) {
    Leaf* const l = (Leaf* const) getRawLeafPointer(n);
    ext_printf("leaf: ");
    for (size_t i = 0; i < l->keyLength; ++i) ext_printf("%c", l->key[i]);
    ext_printf("\n");
    return;
  }
  
  static const char* const sizeNames[] = { "4", "16", "48", "256" };
  
  ext_printf("node %s", sizeNames[n->type]);
  if (n->prefixLength > 0) {
    ext_printf(", partial: '");
    for (size_t i = 0; i < n->prefixLength; ++i) ext_printf("%c", n->partial[i]);
    ext_printf("'");
  }
  ext_printf(", keys: '");
  
#if defined(SUPPRESS_DIAGNOSTIC) && defined(__GNUC__) && !defined(__clang__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
  
  union {
    const Node4  * p1;
    const Node16 * p2;
    const Node48 * p3;
    const Node256* p4;
  } p;
  switch (n->type) {
    case NODE4:
    p.p1 = (const Node4*) n;
    for (uint8_t i = 0; i < n->numChildren; ++i)
      if (p.p1->keys[i] != 0) ext_printf("%c", p.p1->keys[i]);
    break;
    
    case NODE16:
    p.p2 = (const Node16*) n;
    for (uint8_t i = 0; i < n->numChildren; ++i)
      if (p.p2->keys[i] != 0) ext_printf("%c", p.p2->keys[i]);
    break;
    
    case NODE48:
    p.p3 = (const Node48*) n;
    for (uint8_t i = 0; i <= 254; ++i)
      if (p.p3->keys[i] != 0) ext_printf("%c", i);
    break;
    
    case NODE256:
    p.p4 = (const Node256*) n;
    for (uint8_t i = 0; i <= 254; ++i)
      if (p.p4->children[i] != NULL) ext_printf("%c", i);
    break;
  }
  ext_printf("'\n");
  
  switch (n->type) {
    case NODE4:
    for (uint8_t i = 0; i < n->numChildren; ++i) {
      if (p.p1->keys[i] != 0) {
        for (size_t i = 0; i < depth + 1; ++i) ext_printf("  ");
        ext_printf("%c: ", p.p1->keys[i]);
        printAtDepth(p.p1->children[i], depth + 1);
      }
    }
    break;
    
    case NODE16:
    for (uint8_t i = 0; i < n->numChildren; ++i) {
      if (p.p2->keys[i] != 0) {
        for (size_t i = 0; i < depth + 1; ++i) ext_printf("  ");
        ext_printf("%c: ", p.p2->keys[i]);
        printAtDepth(p.p2->children[i], depth + 1);
      }
    }
    break;
    
    case NODE48:
    for (uint8_t i = 0; i <= 254; ++i) {
      if (p.p3->keys[i] != 0) {
        for (size_t j = 0; j < depth + 1; ++j) ext_printf("  ");
        ext_printf("%c: ", i);
        printAtDepth(p.p3->children[p.p3->keys[i] - 1], depth + 1);
      }
    }
    break;
    
    case NODE256:
    for (uint8_t i = 0; i <= 254; ++i) {
      if (p.p4->children[i] != NULL) {
        for (size_t j = 0; j < depth + 1; ++j) ext_printf("  ");
        ext_printf("%c: ", i);
        printAtDepth(p.p4->children[i], depth + 1);
      }
    }
    break;
  }

#if defined(SUPPRESS_DIAGNOSTIC) && defined(__GNUC__) && !defined(__clang__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
#  pragma GCC diagnostic pop
#endif

}

