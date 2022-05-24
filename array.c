#include <assert.h>
#include <string.h>

#include "decs.h"

int array_init(size_t const data_siz, Array ** array) {
    Array * __out = malloc(sizeof(Array));
    RESET(size_t, __out->data_siz, data_siz);
    __out->capacity = 0;
    __out->size     = 0;
    __out->rawdata  = NULL;
    *array = __out;
    return ERROR_NONE;
}

int array_resize(Array * array, unsigned int const newsize,
        void const * val) {
    int ierr;

    if(newsize > array->capacity) {
        ierr = array_reserve(array, newsize); CHKERRQ(ierr);
    }

    unsigned int const oldsize = array->size;
    array->size = newsize;

    if(newsize > oldsize) {
        ierr = array_assign(array, oldsize, newsize, val); CHKERRQ(ierr);
    }

    return ERROR_NONE;
}

int array_reserve(Array * array, unsigned int const capacity) {
    if(capacity > array->capacity) {
        array->rawdata = realloc(array->rawdata, array->data_siz*capacity);
        if(NULL == array->rawdata && capacity > 0) {
            THROW_ERROR(ERROR_OUTOFMEM, "Out of memory");
        }
        array->capacity = capacity;
    }
    return ERROR_NONE;
}

int array_shrink_to_fit(Array * array) {
    if(array->capacity > array->size) {
        array->rawdata = realloc(array->rawdata, array->data_siz*array->size);
        if(NULL == array->rawdata && array->size > 0) {
            THROW_ERROR(ERROR_OUTOFMEM, "Out of memory");
        }
        array->capacity = array->size;
    }
    return ERROR_NONE;
}

int array_at(Array const * array, unsigned int const index,
        void ** elem) {
    assert(index < array->size);
    *elem = array->rawdata + index*array->data_siz;
    return ERROR_NONE;
}

int array_find_first(Array const * array, void const * val,
        void ** elem) {
    char ** p = (char **)elem;
    for(*p = array->rawdata; *p != array->rawdata +
            array->data_siz*array->size; *p += array->data_siz) {
        if(0 == memcmp(val, *p, array->data_siz)) {
            return ERROR_NONE;
        }
    }
    *p = NULL;
    return ERROR_NONE;
}

int array_assign(Array * array, unsigned int const begin,
        unsigned int const end, void const * val) {
    assert(begin <  array->size);
    assert(end   <= array->size);
    for(unsigned int i = begin; i != end; ++i) {
        char * pos = array->rawdata + i*array->data_siz;
        memcpy(pos, val, array->data_siz);
    }
    return ERROR_NONE;
}

int array_clear(Array * array) {
    array->size = 0;
    return ERROR_NONE;
}

int array_push_back(Array * array, void const * val) {
    int ierr;
    if(array->size == array->capacity) {
        unsigned int const newcapacity = MAX(array->capacity + 1,
                (unsigned int)((double)array->capacity*ARRAY_MEMINC));
        ierr = array_reserve(array, newcapacity); CHKERRQ(ierr);
    }
    char * pos = array->rawdata + array->size*array->data_siz;
    memcpy(pos, val, array->data_siz);
    ++array->size;
    return ERROR_NONE;
}

int array_attach(Array * base, Array const * tail) {
    int ierr;
    if(base->data_siz != tail->data_siz) {
        THROW_ERROR(ERROR_LOGICAL, "Trying to concatenate two "
                "arrays of different type");
    }
    if(base->size + tail->size > base->capacity) {
        unsigned int const newcapacity = MAX(base->size + tail->size,
                (unsigned int)((double)base->capacity*ARRAY_MEMINC));
        ierr = array_reserve(base, newcapacity); CHKERRQ(ierr);
    }
    char * pos = base->rawdata + base->size*base->data_siz;
    memcpy(pos, tail->rawdata, tail->data_siz*tail->size);
    base->size = base->size + tail->size;
    return ERROR_NONE;
}

int array_pop_back(Array * array, void * val) {
    --array->size;
    char * pos = array->rawdata + array->size*array->data_siz;
    memcpy(val, pos, array->data_siz);
    return ERROR_NONE;
}

int array_free(Array * array) {
    free(array->rawdata);
    free(array);
    return ERROR_NONE;
}
