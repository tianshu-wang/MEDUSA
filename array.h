//! \file array.h Simple array class

#ifndef ARRAY_H
#define ARRAY_H

#include <stdlib.h>

//! Bare-bone array class
/*!
 *  This class provides a simple implementation of an array of POD objects
 *  stored contiguously in memory.
 *
 *  This class keeps track of its size and of the number of elements stored in
 *  the array. When a new element is added to the array it takes care of
 *  growing the array if needed.
 *
 *  Array works similarly to C++ std::vector. Note, however, that Array
 *  cannot be used to store complex types!
 */
typedef struct {
    //! Size of each of the data elements stored in the array
    size_t const data_siz;
    //! Maximum capacity of the array
    unsigned int capacity;
    //! Number of elements stored in the array
    unsigned int size;
    //! Raw data
    char * rawdata;
} Array;

//! Initialize an array
int array_init(
        //! [in] Size of the objects to store in the array
        size_t const data_siz,
        //! [out] Newly allocated array
        Array ** array);

//! Resize a Array
/*!
 *  This changes the size of the array and, if necessary, also increases
 *  the capacity of the array.
 */
int array_resize(
        //! [in] Array object to resize
        Array * array,
        //! [in] New size of the Array object
        unsigned int const size,
        //! [in] Fill value
        void const * val);

//! Changes the capacity of the array
/*!
 *  This ensures that array->capacity is at least equal to capacity
 */
int array_reserve(
        //! [in,out] Re-sized array
        Array * array,
        //! [in] New size
        unsigned int const capacity);

//! Shrinks the array so that at the end capacity is equal to size
int array_shrink_to_fit(
        //! [in,out] Array to shrink
        Array * array);

//! Gets a reference to an element of the array
int array_at(
        //! [in] Array object
        Array const * array,
        //! [in] Index of the element to access
        unsigned int const index,
        //! [out] Pointer to the wanted element
        void ** elem);

//! Finds the first element having a given value in the array
int array_find_first(
        //! [in] Array object
        Array const * array,
        //! [in] Value of the element to find
        void const * val,
        //! [out] First element found (NULL if none has been found)
        void ** elem);

//! Set all of the elements of the vector to be equal to *val
int array_assign(
        //! [in,out] Array to set
        Array * array,
        //! [in] First element to set
        unsigned const int begin,
        //! [in] Index past the last element to set
        unsigned const int end,
        //! [in] Value
        void const * val);

//! Empties the array
/*!
 *  This sets \e array->size = 0, but does not affect the capacity of the
 *  array
 */
int array_clear(
        //! [in,out] Re-sized array
        Array * array);

//! Appends an element at the end of the array
int array_push_back(
        //! [in,out] Array object
        Array * array,
        //! [in] Element to append to the array
        void const * val);

//! Attach the contents of \e tail at the back of \e base
/*!
 *  This is equivalent to pushing back all of the elements of tail into
 *  the base array, but it is more efficient
 */
int array_attach(
        //! [in,out] Array object
        Array * base,
        //! [in] Array containing the data to append
        Array const * tail);

int array_pop_back(
        //! [in,out] Array object
        Array * array,
        //! [out] Where to store the wanted value
        void * val);

//! De-allocates an Array
int array_free(
        //! [in,out] Array to de-allocate
        Array * array);

//! When growing an array increase the capacity by this factor
#ifndef ARRAY_MEMINC
#define ARRAY_MEMINC      (2)
#endif

#endif

/* vim: set ft=c.doxygen : */
