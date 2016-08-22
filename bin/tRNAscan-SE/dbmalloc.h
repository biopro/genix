/*
 * (c) Copyright 1990, 1991, 1992 Conor P. Cahill (uunet!virtech!cpcahil).  
 * You may copy, distribute, and use this software as long as this
 * copyright statement is not removed.
 */
/*
 * $Id: malloc.h.org,v 1.31 1992/06/22 23:40:10 cpcahil Exp $
 */

#ifndef _DEBUG_MALLOC_INC
#define _DEBUG_MALLOC_INC 1

#ifdef    force_cproto_to_use_defines

/*
 * these are just here because cproto used the c-preprocessor to generate
 * the prototypes and if they were left as #defines the prototypes.h file
 * would have the contents of the define, not the define itself
 */

typedef char		DATATYPE;
typedef int		SIZETYPE;
typedef void		VOIDTYPE;
typedef char		MEMDATA;
typedef int		MEMSIZE;
typedef int		STRSIZE;

#ifdef WRTSIZE
#undef WRTSIZE
#endif
typedef unsigned int	WRTSIZE;

/*
 * for now, define CONST as const.  A sed script in the makefile will change 
 * this back to CONST in the prototypes.h file.
 */
#define CONST const

#else  /* force_cproto_to_use_defines */

/*
 * The following entries are automatically added by the Configure script.
 * If they are not correct for your system, then Configure is not handling
 * your system correctly.  Please report this to the author along with
 * a description of your system and the correct values
 */



#if       __GNUC__ && __STDC__ && __cplusplus

#define VOIDTYPE void
#define CONST const
#define DATATYPE void
#define SIZETYPE size_t
#define MEMDATA void
#define MEMSIZE size_t
#define MEMCMPTYPE unsigned char
#define STRSIZE size_t
#define STRCMPTYPE unsigned char

#else  /* __GNUC__ && __STDC__ && __cplusplus */

#define VOIDTYPE void
#define CONST const
#define DATATYPE void
#define SIZETYPE size_t
#define MEMDATA void
#define MEMSIZE int
#define MEMCMPTYPE unsigned char
#define STRSIZE size_t
#define STRCMPTYPE unsigned char

#endif /* __GNUC__ && __STDC__ && __cplusplus */

/*
 * END of automatic configuration stuff.
 */

/*
 * if DATATYPE is not defined, then the configure script must have had a 
 * problem, or was used with a different compiler.  So we have to stop
 * here and get the user to fix the problem.
 */
#ifndef   DATATYPE
	/*
	 * the following string should cause a comilation error and get the
	 * user to look at this stuff to find out what is wrong.
	 */
	"This file is not configured correctly for this system. Run configure
	 and check its results" 
	char * malloc(); /* DON'T REMOVE THIS LINE if you get a compiler error
			    here it is because the malloc.h file is not 
			    configured correctly  See the readme/problems
			    files for more info */

#endif /* DATATYPE */

#endif /* force_cproto_to_use_defines */

#define VOIDCAST (VOIDTYPE)

/*
 * since we redefine much of the stuff that is #defined in string.h and 
 * memory.h, we should do what we can to make sure that they don't get 
 * included after us.  This is typically accomplished by a special symbol
 * (similar to _DEBUG_MALLOC_INC defined above) that is #defined when the
 * file is included.  Since we don't want the file to be included we will
 * #define the symbol ourselves.  These will typically have to change from
 * one system to another.  I have put in several standard mechanisms used to
 * support this mechanism, so hopefully you won't have to modify this file.
 */
#ifndef _H_STRING
#define _H_STRING		1
#endif 
#ifndef __STRING_H
#define __STRING_H		1
#endif 
#ifndef _STRING_H_
#define _STRING_H_		1
#endif 
#ifndef _STRING_H 
#define _STRING_H 		1
#endif 
#ifndef _STRING_INCLUDED
#define _STRING_INCLUDED	1
#endif
#ifndef __string_h
#define __string_h		1
#endif
#ifndef _string_h
#define _string_h		1
#endif
#ifndef _strings_h
#define _strings_h		1
#endif
#ifndef __strings_h
#define __strings_h		1
#endif
#ifndef _H_MEMORY
#define _H_MEMORY		1
#endif
#ifndef __MEMORY_H
#define __MEMORY_H		1
#endif
#ifndef _MEMORY_H_
#define _MEMORY_H_		1
#endif
#ifndef _MEMORY_H
#define _MEMORY_H		1
#endif
#ifndef _MEMORY_INCLUDED
#define _MEMORY_INCLUDED	1
#endif
#ifndef __memory_h
#define __memory_h		1
#endif
#ifndef _memory_h
#define _memory_h		1
#endif

/*
 * Malloc warning/fatal error handler defines...
 */
#define M_HANDLE_DUMP	0x80  /* 128 */
#define M_HANDLE_IGNORE	0
#define M_HANDLE_ABORT	1
#define M_HANDLE_EXIT	2
#define M_HANDLE_CORE	3
	
/*
 * Mallopt commands and defaults
 *
 * the first four settings are ignored by the debugging mallopt, but are
 * here to maintain compatibility with the system malloc.h.
 */
#define M_MXFAST	1		/* ignored by mallopt		*/
#define M_NLBLKS	2		/* ignored by mallopt		*/
#define M_GRAIN		3		/* ignored by mallopt		*/
#define M_KEEP		4		/* ignored by mallopt		*/
#define MALLOC_WARN	100		/* set malloc warning handling	*/
#define MALLOC_FATAL	101		/* set malloc fatal handling	*/
#define MALLOC_ERRFILE	102		/* specify malloc error file	*/
#define MALLOC_CKCHAIN	103		/* turn on chain checking	*/
#define MALLOC_FILLAREA	104		/* turn off area filling	*/
#define MALLOC_LOWFRAG	105		/* use best fit allocation mech	*/
#define MALLOC_CKDATA	106		/* turn off/on data checking	*/
#define MALLOC_REUSE	107		/* turn off/on freed seg reuse	*/
#define MALLOC_SHOWLINKS 108		/* turn off/on adjacent link disp */

union malloptarg
{
	int	  i;
	char	* str;
};

/*
 * Malloc warning/fatal error codes
 */

#define M_CODE_CHAIN_BROKE	1	/* malloc chain is broken	*/
#define M_CODE_NO_END		2	/* chain end != endptr		*/
#define M_CODE_BAD_PTR		3	/* pointer not in malloc area	*/
#define M_CODE_BAD_MAGIC	4	/* bad magic number in header	*/
#define M_CODE_BAD_CONNECT	5	/* chain poingers corrupt	*/
#define M_CODE_OVERRUN		6	/* data overrun in malloc seg	*/
#define M_CODE_REUSE		7	/* reuse of freed area		*/
#define M_CODE_NOT_INUSE	8	/* pointer is not in use	*/
#define M_CODE_NOMORE_MEM	9	/* no more memory available	*/
#define M_CODE_OUTOF_BOUNDS	10	/* gone beyound bounds 		*/
#define M_CODE_FREELIST_BAD	11	/* inuse segment on freelist	*/
#define M_CODE_NOBOUND		12	/* can't calculate boundry	*/
#define M_CODE_STK_NOCUR	13	/* no current element on stack	*/
#define M_CODE_STK_BADFUNC	14	/* current func doesn't match	*/
#define M_CODE_UNDERRUN		15	/* data underrun in malloc seg	*/

#ifndef __STDCARGS
#if  __STDC__ || __cplusplus
#define __STDCARGS(a) a
#else
#define __STDCARGS(a) ()
#endif
#endif

#if __cplusplus
extern "C" {
#endif

VOIDTYPE	  malloc_dump __STDCARGS((int));
VOIDTYPE	  malloc_list __STDCARGS((int,unsigned long, unsigned long));
int		  mallopt __STDCARGS((int, union malloptarg));
DATATYPE	* debug_calloc __STDCARGS((CONST char *,int,SIZETYPE,SIZETYPE));
VOIDTYPE	  debug_cfree __STDCARGS((CONST char *, int, DATATYPE *));
VOIDTYPE	  debug_free __STDCARGS((CONST char *, int, DATATYPE *));
DATATYPE	* debug_malloc __STDCARGS((CONST char *,int, SIZETYPE));
DATATYPE	* debug_realloc __STDCARGS((CONST char *,int,
					    DATATYPE *,SIZETYPE));
unsigned long	  DBmalloc_size __STDCARGS((CONST char *,int,unsigned long *));
int		  DBmalloc_chain_check __STDCARGS((CONST char *,int,int));

void		  StackEnter __STDCARGS((CONST char *, CONST char *, int));
void		  StackLeave __STDCARGS((CONST char *, CONST char *, int));

/*
 * X allocation related prototypes
 */
char		* debug_XtMalloc __STDCARGS((CONST char *, int, unsigned int));
char		* debug_XtRealloc __STDCARGS((CONST char *, int,
						char *, unsigned int));
char		* debug_XtCalloc __STDCARGS((CONST char *, int,
						unsigned int, unsigned int));
void		  debug_XtFree __STDCARGS((CONST char *, int, char *));
void		* debug_XtBCopy  __STDCARGS((CONST char *, int, char *,
						char *, int));
extern void	(*XtAllocErrorHandler) __STDCARGS((CONST char *));

/*
 * memory(3) related prototypes
 */
MEMDATA  	* DBmemccpy __STDCARGS((CONST char *file, int line,
					MEMDATA  *ptr1, CONST MEMDATA  *ptr2,
					int ch, MEMSIZE len));
MEMDATA  	* DBmemchr __STDCARGS((CONST char *file, int line,
					CONST MEMDATA  *ptr1, int ch,
					MEMSIZE len));
MEMDATA 	* DBmemmove __STDCARGS((CONST char *file, int line,
					MEMDATA  *ptr1, CONST MEMDATA  *ptr2,
					MEMSIZE len));
MEMDATA 	* DBmemcpy __STDCARGS((CONST char *file, int line,
					MEMDATA  *ptr1, CONST MEMDATA  *ptr2,
					MEMSIZE len));
int		  DBmemcmp __STDCARGS((CONST char *file, int line,
					CONST MEMDATA  *ptr1,
					CONST MEMDATA  *ptr2, MEMSIZE len));
MEMDATA 	* DBmemset __STDCARGS((CONST char *file, int line,
					MEMDATA  *ptr1, int ch, MEMSIZE len));
MEMDATA 	* DBbcopy __STDCARGS((CONST char *file, int line,
					CONST MEMDATA  *ptr2, MEMDATA  *ptr1,
					MEMSIZE len));
MEMDATA  	* DBbzero __STDCARGS((CONST char *file, int line,
					MEMDATA  *ptr1, MEMSIZE len));
int		  DBbcmp __STDCARGS((CONST char *file, int line,
					CONST MEMDATA  *ptr2,
					CONST MEMDATA  *ptr1, MEMSIZE len));

/*
 * string(3) related prototypes
 */
char		* DBstrcat __STDCARGS((CONST char *file,int line, char *str1,
					CONST char *str2));
char		* DBstrdup __STDCARGS((CONST char *file, int line,
					CONST char *str1));
char		* DBstrncat __STDCARGS((CONST char *file, int line, char *str1,
					CONST char *str2, STRSIZE len));
int		  DBstrcmp __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2));
int		  DBstrncmp __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2,
					STRSIZE len));
int		  DBstricmp __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2));
int		  DBstrincmp __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2,
					STRSIZE len));
char		* DBstrcpy __STDCARGS((CONST char *file, int line, char *str1,
					CONST char *str2));
char		* DBstrncpy __STDCARGS((CONST char *file, int line, char *str1,
					CONST char *str2, STRSIZE len));
STRSIZE		  DBstrlen __STDCARGS((CONST char *file, int line,
					CONST char *str1));
char		* DBstrchr __STDCARGS((CONST char *file, int line,
					CONST char *str1, int c));
char		* DBstrrchr __STDCARGS((CONST char *file, int line,
					CONST char *str1, int c));
char		* DBindex __STDCARGS((CONST char *file, int line,
					CONST char *str1, int c));
char		* DBrindex __STDCARGS((CONST char *file, int line,
					CONST char *str1, int c));
char		* DBstrpbrk __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2));
STRSIZE		  DBstrspn __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2));
STRSIZE		  DBstrcspn __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2));
char		* DBstrstr __STDCARGS((CONST char *file, int line,
					CONST char *str1, CONST char *str2));
char		* DBstrtok __STDCARGS((CONST char *file, int line, char *str1,
					CONST char *str2));

#if __cplusplus
};
#endif

/*
 * Macro which enables logging of the file and line number for each allocation
 * so that it is easier to determine where the offending malloc comes from.
 *
 * NOTE that only code re-compiled with this include file will have this 
 * additional info.  Calls from libraries that have not been recompiled will
 * just have a null string for this info.
 */
#ifndef IN_MALLOC_CODE

/*
 * allocation functions
 */
#define malloc(len)		debug_malloc( __FILE__,__LINE__, (len))
#define realloc(ptr,len)	debug_realloc(__FILE__,__LINE__, (ptr), (len))
#define calloc(numelem,size)	debug_calloc(__FILE__,__LINE__,(numelem),(size))
#define cfree(ptr)		debug_cfree(__FILE__,__LINE__,(ptr))
#define free(ptr)		debug_free(__FILE__,__LINE__,(ptr))
#define malloc_size(histptr)	DBmalloc_size(__FILE__,__LINE__,(histptr))
#define malloc_chain_check(todo) DBmalloc_chain_check(__FILE__,__LINE__,(todo))

/* 
 * X allocation routines
 */
#define XtCalloc(_num,_size)	debug_XtCalloc(__FILE__,__LINE__,_num,_size)
#define XtMalloc(_size)		debug_XtMalloc(__FILE__,__LINE__,_size)
#define XtRealloc(_ptr,_size)	debug_XtRealloc(__FILE__,__LINE__,_ptr,_size)
#define XtFree(_ptr)		debug_XtFree(__FILE__,__LINE__,_ptr)
#define _XtBCopy(ptr1,ptr2,len) debug_XtBcopy(__FILE__,__LINE__,ptr1,ptr2,len)

/*
 * Other allocation functions
 */
#define _malloc(_size)		debug_malloc(__FILE__,__LINE__,_size)
#define _realloc(_ptr,_size)	debug_realloc(__FILE__,__LINE__,_ptr,_size)
#define _calloc(_num,_size)	debug_calloc(__FILE__,__LINE__,_num,_size)
#define _free(_ptr)		debug_free(__FILE__,__LINE__,_ptr)

/*
 * memory(3) related functions
 */
#ifdef bcopy
#undef bcopy
#endif
#ifdef bzero
#undef bzero
#endif
#ifdef bcmp
#undef bcmp
#endif
#define memccpy(ptr1,ptr2,ch,len) DBmemccpy(__FILE__,__LINE__,ptr1,ptr2,ch,len)
#define memchr(ptr1,ch,len)	  DBmemchr(__FILE__,__LINE__,ptr1,ch,len)
#define memmove(ptr1,ptr2,len)    DBmemmove(__FILE__,__LINE__,ptr1, ptr2, len)
#define memcpy(ptr1,ptr2,len)     DBmemcpy(__FILE__, __LINE__, ptr1, ptr2, len)
#define memcmp(ptr1,ptr2,len)     DBmemcmp(__FILE__,__LINE__,ptr1, ptr2, len)
#define memset(ptr1,ch,len)       DBmemset(__FILE__,__LINE__,ptr1, ch, len)
#define bcopy(ptr2,ptr1,len)      DBbcopy(__FILE__,__LINE__,ptr2,ptr1,len)
#define bzero(ptr1,len)           DBbzero(__FILE__,__LINE__,ptr1,len)
#define bcmp(ptr2,ptr1,len)       DBbcmp(__FILE__, __LINE__, ptr2, ptr1, len)

#define _bcopy(ptr2,ptr1,len)     DBbcopy(__FILE__,__LINE__,ptr2,ptr1,len)
#define _bzero(ptr1,len)          DBbzero(__FILE__,__LINE__,ptr1,len)
#define _bcmp(ptr2,ptr1,len)      DBbcmp(__FILE__,__LINE__,ptr2,ptr1,len)
#define __dg_bcopy(ptr2,ptr1,len) DBbcopy(__FILE__,__LINE__,ptr2,ptr1,len)
#define __dg_bzero(ptr1,len)      DBbzero(__FILE__,__LINE__,ptr1,len)
#define __dg_bcmp(ptr2,ptr1,len)  DBbcmp(__FILE__,__LINE__,ptr2,ptr1,len)

/*
 * string(3) related functions
 */
#ifdef index
#undef index
#endif
#ifdef rindex
#undef rindex
#endif
#ifdef strcpy
#undef strcpy
#endif
#ifdef strcpy
#undef strcmp
#endif
#define index(str1,c)		  DBindex(__FILE__, __LINE__, str1, c)
#define rindex(str1,c)		  DBrindex(__FILE__, __LINE__, str1, c)
#define strcat(str1,str2)	  DBstrcat(__FILE__,__LINE__,str1,str2)
#define strchr(str1,c)		  DBstrchr(__FILE__, __LINE__, str1,c)
#define strcmp(str1,str2)	  DBstrcmp(__FILE__, __LINE__, str1, str2)
#define strcpy(str1,str2)	  DBstrcpy(__FILE__, __LINE__, str1, str2)
#define strcspn(str1,str2)	  DBstrcspn(__FILE__, __LINE__, str1,str2)
#define strdup(str1)		  DBstrdup(__FILE__, __LINE__, str1)
#define stricmp(str1,str2)	  DBstricmp(__FILE__, __LINE__, str1, str2)
#define strincmp(str1,str2,len)	  DBstrincmp(__FILE__, __LINE__, str1,str2,len)
#define strlen(str1)		  DBstrlen(__FILE__, __LINE__, str1)
#define strncat(str1,str2,len)	  DBstrncat(__FILE__, __LINE__, str1,str2,len)
#define strncpy(str1,str2,len)	  DBstrncpy(__FILE__,__LINE__,str1,str2,len)
#define strncmp(str1,str2,len)	  DBstrncmp(__FILE__, __LINE__, str1,str2,len)
#define strpbrk(str1,str2)	  DBstrpbrk(__FILE__, __LINE__, str1,str2)
#define strrchr(str1,c)		  DBstrrchr(__FILE__,__LINE__,str1,c)
#define strspn(str1,str2)	  DBstrspn(__FILE__, __LINE__, str1,str2)
#define strstr(str1,str2)	  DBstrstr(__FILE__, __LINE__, str1, str2)
#define strtok(str1,str2)	  DBstrtok(__FILE__, __LINE__, str1, str2)

/*
 * malloc stack related functions
 */
#define malloc_enter(func)	  StackEnter(func,__FILE__,__LINE__)
#define malloc_leave(func)	  StackLeave(func,__FILE__,__LINE__)

#endif /* IN_MALLOC_CODE */

#endif /* _DEBUG_MALLOC_INC */

/*
 * $Log: malloc.h.org,v $
 * Revision 1.31  1992/06/22  23:40:10  cpcahil
 * many fixes for working on small int systems
 *
 * Revision 1.30  1992/05/06  04:53:29  cpcahil
 * performance enhancments
 *
 * Revision 1.29  1992/04/22  18:17:32  cpcahil
 * added support for Xt Alloc functions, linted code
 *
 * Revision 1.28  1992/04/13  19:08:18  cpcahil
 * fixed case insensitive stuff
 *
 * Revision 1.27  1992/04/13  18:41:18  cpcahil
 * added case insensitive string comparison routines
 *
 * Revision 1.26  1992/04/13  17:26:25  cpcahil
 * minor portability changes
 *
 * Revision 1.25  1992/04/13  14:13:18  cpcahil
 * cleanup of log message.
 *
 * Revision 1.24  1992/04/13  03:09:14  cpcahil
 * lots of changes.
 *
 * Revision 1.23  1992/03/01  12:42:38  cpcahil
 * added support for managing freed areas and fixed doublword bndr problems
 *
 * Revision 1.22  1992/02/07  15:51:07  cpcahil
 * mods for sun4
 *
 * Revision 1.21  1992/01/29  01:35:32  cpcahil
 * added sgi definition.
 *
 * Revision 1.20  1992/01/28  21:42:25  cpcahil
 * changes for the ibmRS6000
 *
 * Revision 1.19  1992/01/28  18:05:37  cpcahil
 * misc fixes for patch 7
 *
 * Revision 1.18  1992/01/22  16:21:35  cpcahil
 * added code to prevent inclusions of string.h and memory.h after malloc.h
 * was included.
 *
 * Revision 1.17  1992/01/10  17:26:46  cpcahil
 * fixed prototypes use of void.
 *
 * Revision 1.16  1992/01/10  16:53:39  cpcahil
 * added more info on sizetype and datatype. added support for overriding
 * use of void type.
 *
 * Revision 1.15  1992/01/09  17:19:11  cpcahil
 * put the close brace in the correct position.
 *
 * Revision 1.14  1992/01/09  17:12:36  cpcahil
 * added code to support inclusion in C++ modules
 *
 * Revision 1.13  1991/12/31  21:31:26  cpcahil
 * changes for patch 6.  See CHANGES file for more info
 *
 * Revision 1.12  1991/12/26  22:31:29  cpcahil
 * added check to make sure file is not included twice.
 *
 * Revision 1.11  1991/12/06  17:58:46  cpcahil
 * added cfree() for compatibility with some wierd systems
 *
 * Revision 1.10  91/12/06  08:54:18  cpcahil
 * cleanup of __STDC__ usage and addition of CHANGES file
 * 
 * Revision 1.9  91/12/04  09:23:40  cpcahil
 * several performance enhancements including addition of free list
 * 
 * Revision 1.8  91/12/02  19:10:11  cpcahil
 * changes for patch release 5
 * 
 * Revision 1.7  91/11/25  14:42:00  cpcahil
 * Final changes in preparation for patch 4 release
 * 
 * Revision 1.6  91/11/24  00:49:28  cpcahil
 * first cut at patch 4
 * 
 * Revision 1.5  91/11/20  11:54:10  cpcahil
 * interim checkin
 * 
 * Revision 1.4  90/08/29  22:23:38  cpcahil
 * fixed mallopt to use a union as an argument.
 * 
 * Revision 1.3  90/05/11  11:04:10  cpcahil
 * took out some extraneous lines
 * 
 * Revision 1.2  90/05/11  00:13:09  cpcahil
 * added copyright statment
 * 
 * Revision 1.1  90/02/23  07:09:03  cpcahil
 * Initial revision
 * 
 */
