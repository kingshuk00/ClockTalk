/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */

/**
 * @file paraver_file.h
 * @brief General-purpose header-only library to aid read Paraver files
 *
 * A Paraver trace file has three sections in order:dfd
 *
 *   1. **Header line**: a single line beginning with `#Paraver (` that encodes
 *      the trace duration, time unit, node count, application count, process
 *      counts per application, and communicator counts.
 *
 *   2. **Communicator section**: one line per communicator, each of the form
 *      `c:<app>:<id>:<size>:<rank0>:<rank1>:...`, listing the global MPI ranks
 *      that belong to each communicator.
 *
 *   3. **Records section**: the body of the trace.  Each line begins with a
 *      record-type digit: `1` = state, `2` = event, `3` = point-to-point
 *      message.
 *
 * The consumer of this library calls `PrvFile_process()` to process every line
 * of the trace using its own callback.
 *   - This function requires the callback to be set before it is called.
 *
 * The consumer can implement a multi-pass workflow by calling the
 * `PrvFile_reloadRecords()` function.
 *   - This function seeks back to the start of the records sections, from where
 *     `PrvFile_process()` can be called once again.
 *
 * The records section is read in 32 MB chunks using `fread()` for throughput.
 * Partial lines at chunk boundaries are carried over to the next chunk via
 * `memmove()`. The time spent in `fread()` is measured separately from processing
 * time and returned by the `PrvFile_process()` function to the consumer.
 *
 * @note Requires `_LARGEFILE_SOURCE` for `fseeko()`/`ftello()` on 32-bit
 *       systems so that file offsets are 64-bit and traces larger than 2 GB
 *       are handled correctly. On 64-bit systems this has no effect.
 *
 * @note Multi-application traces are not fully supported. A warning is
 *       printed if numApps > 1 and results may be inaccurate.
 *
 * @note Designed to be used standalone.
 */

#ifndef CLOCKTALK_PARAVER_PARAVER_FILE_H__
#define CLOCKTALK_PARAVER_PARAVER_FILE_H__

#define _LARGEFILE_SOURCE

/******************************************************************************/
/* Public APIs of this header-only library                                    */
/******************************************************************************/

typedef struct ParaverFile_struct__ ParaverFile;

/**
 * @brief Opens a Paraver trace file and parses its header.
 *
 * Allocates and returns a `ParaverFile` handle populated with all metadata
 * from the header and communicator sections. The file pointer is left
 * positioned at the start of the records section.
 *
 * Specifically:
 *   - Opens the file and records its total size via `fstat()`.
 *   - Reads and parses the single header line to extract runtime, time unit,
 *     node count, application count, process counts, and communicator counts.
 *   - Reads each communicator line to accumulate the total membership count
 *     in @p `allCommsSizes` (sizes are not stored yet; use `PrvFile_readComms()`
 *     for that).
 *   - Records file offsets for the communicator section and records section so
 *     they can be seeked to independently.
 *
 * @param[in] fn Path to the `.prv` trace file.
 * @return Allocated and initialised ParaverFile handle, or `NULL` on any error
 *         (eg. file not found, unrecognised header format, allocation failure).
 *          All resources are cleaned up before returning `NULL`.
 */
inline static ParaverFile *PrvFile_open(const char *const fn);

/**
 * @brief Closes the file handle and frees the `ParaverFile` struct.
 *
 * Safe to call with `NULL`.
 *
 * @param[in] file Handle to close. Must not be used after this call.
 */
inline static void PrvFile_close(ParaverFile *const file);

/**
 * @brief Returns the trace duration in the trace's native time unit.
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return Trace duration as parsed from the header.
 */
inline static long long PrvFile_runTime(const ParaverFile *const file);

/**
 * @brief Returns the time unit string declared in the trace header.
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return "ns" if the header declares nanoseconds, "us" otherwise.
 */
inline static const char *PrvFile_timeUnit(const ParaverFile *const file);

/**
 * @brief Returns the number of hardware nodes recorded in the header.
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return Node count.
 */
inline static int PrvFile_numNodes(const ParaverFile *const file);

/**
 * @brief Returns the number of applications recorded in the header.
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return Application count. Values > 1 indicate a multi-app trace, which
 *         is not fully supported.
 */
inline static int PrvFile_numApps(const ParaverFile *const file);

/**
 * @brief Returns the total number of MPI processes across all applications.
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return Process count.
 */
inline static int PrvFile_numProcs(const ParaverFile *const file);

/**
 * @brief Returns the number of communicators declared in the trace.
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return Communicator count.
 */
inline static int PrvFile_numComms(const ParaverFile *const file);

/**
 * @brief Returns the total number of communicator membership entries.
 *
 * This is the sum of sizes of all communicators and the number of entries
 * needed in the flat ranks array when allocating communicator information.
 *
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return Total membership entry count across all communicators.
 */
inline static int PrvFile_allCommsSizes(const ParaverFile *const file);

/**
 * @brief Reads communicator membership data into caller-supplied arrays.
 *
 * Seeks to the communicator section and reads each communicator line,
 * populating @p commsSizes with the size of each communicator and
 * @p commsRanks with the global MPI ranks of its members.
 *
 * The caller is responsible for allocating the arrays. The expected layout
 * is a flat contiguous block for all ranks with @p commsRanks[i] pointing
 * into it: `commsRanks[0]` = flat_block; `commsRanks[i]` = `commsRanks[i-1] +
 * commsSizes[i-1]`.
 *
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @param[inout] commsSizes Array of length numComms. `commsSizes[i]` will be set
 *                          to the number of members of i-th communicator.
 * @param[inout] commsRanks Array of pointers of length numComms.
 *                          `commsRanks[i]` points to the start of the rank list
 *                          for communicator i. Ranks are 0-based.
 * @return 0 on success, -1 on seek failure.
 */
inline static int PrvFile_readComms(const ParaverFile *const file,
                                   int *const commsSizes,
                                   int **const commsRanks);

/**
 * @brief Seeks the file pointer back to the start of the records section.
 *
 * Possible use is between a count pass and a read pass at the consumer to
 * restart record processing without reopening the file.
 *
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @return 0 on success, -1 on seek failure.
 */
inline static int PrvFile_reloadRecords(const ParaverFile *const file);

/**
 * @brief Sets the line callback invoked for each record line during processing.
 *
 * The callback receives a null-terminated string containing one record line
 * with the trailing newline replaced by '\0'. It is called once per line
 * during `PrvFile_process()`.
 *
 * Must be called before `PrvFile_process()` to set the callback.
 * Call it again with a different function between passes.
 *
 * @param[in] file Opened file handle returned by `PrvFile_open()`.
 * @param[in] processorFunc Callback function that processes one mutable record
 *                          line at a time. May modify the string in place
 *                          (e.g. via `strtok()`) but must not free it.
 */
inline static void PrvFile_setLineProcessor(ParaverFile *const file,
                                            void (*processorFunc)(char *const));

/**
 * @brief Reads and processes the entire records section.
 *
 * Reads the file in 32 MB chunks using `fread()` for throughput. Partial lines
 * at chunk boundaries are carried forward with `memmove()`. The line processor
 * set by `PrvFile_setLineProcessor()` is called once for each complete line.
 *
 * The time spent inside fread() is measured separately from processing time
 * using `CLOCK_MONOTONIC` and returned to the caller for I/O performance
 * reporting.
 *
 * @param[in] file Opened file handle returned by `PrvFile_open()`. Must have a
                   line processor set via `PrvFile_setLineProcessor()`.
 * @param[in] verbosity If non-zero, prints a progress percentage to `stdout` as
 *                      processing proceeds.
 * @return Time spent in `fread()` in seconds. Returns 0.0 if no line processor
 *         has been set.
 */
inline static double PrvFile_process(ParaverFile *const file,
                                     const int verbosity);

/**
 * @brief Advances a record line pointer past the next ':' field separator.
 *
 * Paraver record lines use ':' as the field delimiter. This function finds
 * the next ':' in the string and returns a pointer to the character
 * immediately after it, i.e. to the start of the next field value.
 *
 * @param[in] p Pointer into a record line, positioned at or before a ':'.
 * @return Pointer to the first character of the next field.
 * @note Asserts that a ':' is found.  Passing a pointer past the last
 *       field will trigger an assertion failure.
 */
inline static char *PrvFile_nextRecNum(char *const p);

/**
 * @brief Advances a record line pointer past @p `n` ':' field separators.
 *
 * Convenience wrapper around `PrvFile_nextRecNum()` for skipping multiple
 * fields at once.
 *
 * @param[in] p Pointer into a record line.
 * @param[in] n Number of ':' separators to skip.
 * @return Pointer to the first character of the field after the nth
 *         separator.
 */
inline static char *PrvFile_nthRecNum(char *const p, const int n);

/**
 * @brief Returns the human-readable name of an MPI event by its Paraver id.
 *
 * Paraver encodes MPI functions as integer event values in the range [0, 216].
 * This function maps those values to their standard MPI function name strings
 * (e.g. 1: "Send", 7: "Bcast"). Event id 0 is "Useful" (out-of-MPI).
 *
 * @param[in] eventId Paraver MPI event id.
 * @return Pointer to a string literal with the function name, or NULL if
 *         @p `eventId` is out of the valid range.
 */
inline static const char *PrvFile_MPIName(const int eventId);

/* -------------------------------------------------------------------------- */
/* END of public APIs of this library                                         */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/* Private APIs of this library                                               */
/* -------------------------------------------------------------------------- */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include<time.h>
#include<limits.h>
#include<assert.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>

/**
 * @brief Handle for an opened Paraver trace file.
 *
 * Populated entirely by `PrvFile_open()`. All fields are read-only after
 * construction; only the lineProcessor is mutated during use.
 *
 * @var ParaverFile::fp
 *   Opened file handle. Positioned at the start of the records section after
 *   `PrvFile_open()` returns, then repositioned by `PrvFile_reloadRecords()`.
 *
 * @var ParaverFile::size
 *   Total file size in bytes, obtained via `fstat()`. Used by
 *   `PrvFile_process()` to compute and display progress as a percentage.
 *
 * @var ParaverFile::commsPos
 *   File offset of the first communicator line, i.e. the byte immediately
 *   after the header line. Used by `PrvFile_readComms()` to seek back to
 *   the communicator section.
 *
 * @var ParaverFile::recsPos
 *   File offset of the first record line, i.e. the byte immediately after the
 *   last communicator line. Used by PrvFile_reloadRecords() for multi-pass
 *   processing.
 *
 * @var ParaverFile::runTime
 *   Trace duration in the time unit given by @p `timeUnit`, as parsed from the
 *   header.
 *
 * @var ParaverFile::timeUnit
 *   Null-terminated string: "ns" if the header contains "_ns", "us" otherwise.
 *
 * @var ParaverFile::numNodes
 *   Number of hardware nodes recorded in the header.
 *
 * @var ParaverFile::numApps
 *   Number of applications recorded in the header. Values > 1 are not fully
 *   supported.
 *
 * @var ParaverFile::numProcs
 *   Total number of MPI processes across all applications.
 *
 * @var ParaverFile::numComms
 *   Number of communicators declared in the header and communicator section.
 *
 * @var ParaverFile::allCommsSizes
 *   Sum of all communicator sizes, i.e. the total number of
 *   (communicator: ranks) membership entries. Used to allocate the flat ranks
 *   array.
 *
 * @var ParaverFile::lineProcessor
 *   Callback set by `PrvFile_setLineProcessor()`, called once per record line
 *   during `ParaverFileProcess()`.
 */
typedef struct ParaverFile_struct__ {
  FILE *fp;
  off_t size;
  off_t commsPos;
  off_t recsPos;
  int verbosity;

  int64_t runTime;
  char timeUnit[4];
  int numNodes;
  int numApps;
  int numProcs;
  int numComms;
  long allCommsSizes;

  void (*lineProcessor)(char *const);
} ParaverFile;


inline static long long PrvFile_runTime(const ParaverFile *const file)
{
  return file->runTime;
}
inline static const char *PrvFile_timeUnit(const ParaverFile *const file)
{
  return file->timeUnit;
}
inline static int PrvFile_numNodes(const ParaverFile *const file)
{
  return file->numNodes;
}
inline static int PrvFile_numApps(const ParaverFile *const file)
{
  return file->numApps;
}
inline static int PrvFile_numProcs(const ParaverFile *const file)
{
  return file->numProcs;
}
inline static int PrvFile_numComms(const ParaverFile *const file)
{
  return file->numComms;
}
inline static int PrvFile_allCommsSizes(const ParaverFile *const file)
{
  return file->allCommsSizes;
}

inline static void PrvFile_setLineProcessor(ParaverFile *const file,
                                            void (*processorFunc)(char *const))
{
  file->lineProcessor= processorFunc;
}

inline static char *PrvFile_nextRecNum(char *const p)
{
  char *x= strchr(p, ':');
  assert(NULL!= x);
  return x+ 1;
}

inline static char *PrvFile_nthRecNum(char *const p, const int n)
{
  char *ptr= p;
  for(int i= 0; i< n; ++i) {
    ptr= PrvFile_nextRecNum(ptr);
  }
  return ptr;
}

inline static off_t prvfile_getSize(FILE *fp)
{
  struct stat fpStat;
  if(fstat(fileno(fp), &fpStat)< 0) {
    fprintf(stderr, "%s: Error reading file statistics\n", __func__);
    return -1;
  }

  if(-1== fseeko(fp, 0, SEEK_SET)) {
    fprintf(stderr, "%s: Error locating Paraver header section.\n", __func__);
    return -1;
  }

  return fpStat.st_size;
}

inline static off_t prvfile_readHeader(FILE *fp, char **headerp, size_t *headerLenp)
{
  ssize_t linelen= getline(headerp, headerLenp, fp);
  if(-1== linelen) {
    if(NULL!= (*headerp)) {
      free(*headerp);
      *headerp= NULL;
    }
    return -1;
  } else {
    *headerLenp= (size_t) linelen;
  }
  return ftello(fp);
}
inline static ParaverFile *PrvFile_open(const char *const fn)
{
  FILE *fp= NULL;
  char *header= NULL;
  ParaverFile *file= (ParaverFile *) malloc(sizeof(ParaverFile));
  if(NULL== file) {
    fprintf(stderr, "%s: Error allocating memroy\n", __func__);
    goto bad;
  }
  memset(file, 0, sizeof(ParaverFile));

  fp= fopen(fn, "r");
  if(NULL== fp) {
    fprintf(stderr, "%s: Error opening file-\"%s\"\n", __func__, fn);
    goto bad;
  }

  if((file->size= prvfile_getSize(fp))< 0) {
    goto bad;
  }

  size_t headerLen= 0;
  if((file->commsPos= prvfile_readHeader(fp, &header, &headerLen))< 0) {
    goto bad;
  }

  char *ptr= header;
  if(0!= strncmp("#Paraver (", ptr, 10)) {
    fprintf(stderr, "%s: Unexpected Paraver trace format - invalid header.\n", __func__);
  }

  ptr= strchr(ptr, ')')+ 2;
  file->runTime= (int64_t) atoll(ptr);

  if(0== strncmp("_ns", ptr- 3, 3)) {
    strcpy(file->timeUnit, "ns");
  } else {
    strcpy(file->timeUnit, "us");
  }

  ptr= PrvFile_nextRecNum(ptr);
  file->numNodes= atoi(ptr);

  ptr= PrvFile_nextRecNum(ptr);
  file->numApps= atoi(ptr);
  if(file->numApps> 1) {
    fprintf(stderr, "Paraver traces with more than 1 applications is not yet "
            "supported.\nResults will be inaccurate.\n");
  }

  file->numProcs= 0;
  file->numComms= 0;
  ptr= PrvFile_nextRecNum(ptr);
  while(NULL!= ptr) {
    file->numProcs+= atoi(ptr);
    ptr= strchr(ptr, '(')+ 1;
    ptr= strchr(ptr, ')')+ 1;
    if(strlen(ptr)> 1) {
      ++ptr;
      file->numComms+= atoi(ptr);
    }
    if(NULL!= strchr(ptr, ':')) {
      ptr= strchr(ptr, ':');
    } else if(NULL!= strchr(ptr, ',')) {
      ptr= strchr(ptr, ',');
    } else {
      ptr= NULL;
    }
  }

  for(int i= 0; i< file->numComms; ++i) {
    ssize_t len= getline(&header, &headerLen, fp);
    if(-1== len) {
      fprintf(stderr, "%s: getline() failed.\n", __func__);
      goto bad;
    }

    int cnp;
    /*                        c:app: id:np:p0:p1... */
    if(1!= sscanf(header, "c:%*d:%*d:%d:", &cnp)) {
      fprintf(stderr, "%s: Unexpected Paraver trace format - invalid communicators.\n",
              __func__);
    }

    (file->allCommsSizes)+= cnp;
  }
  file->recsPos= ftello(fp);

  file->fp= fp;
  goto bye;

bad:
  if(NULL!= fp) {
    fclose(fp);
    fp= NULL;
  }
  if(NULL!= file) {
    free(file);
    file= NULL;
  }

bye:
  if(NULL!= header) {
    free(header);
    header= NULL;
  }
  return file;
}

inline static void PrvFile_close(ParaverFile *const file)
{
  if(NULL!= file) {
    if(NULL!= file->fp) {
      fclose(file->fp);
      file->fp= NULL;
    }
    free(file);
  }
}

inline static int prvFile_readOneComm(char *const str, const int ix,
                                      int *const cs, int *const crs)
{
  int cix;
  /*                ignored-c:app: */
  if(2!= sscanf(str, "c:%*d:%d:%d:", &cix, cs)) {
    fprintf(stderr,
            "%s: Unexpected Paraver trace format - invalid communicator.\n%s\n",
            __func__, str);
  }
  --cix;
  if(cix!= ix) {
    fprintf(stderr, "%s: Erratic communicator index\n%s\n", __func__, str);
  }
  char *ptr= PrvFile_nthRecNum(str, 3);
  for(int i= 0; i< *cs; ++i) {
    ptr= PrvFile_nextRecNum(ptr);
    crs[i]= atoi(ptr)- 1;
  }
  return *cs;
}

inline static int PrvFile_readComms(const ParaverFile *const file,
                                    int *const cs, int **const crs)
{
  if(file->numComms< 1) {
    return 0;
  }

  if(-1== fseeko(file->fp, file->commsPos, SEEK_SET)) {
    fprintf(stderr, "%s: Cannot rewind communicators in file!\n", __func__);
    return -1;
  }

  char *str; size_t n= 0;
  ssize_t len= getline(&str, &n, file->fp);
  if(-1== len) {
    fprintf(stderr, "%s: Call to getline() failed.\n", __func__);
  }
  int ncps= prvFile_readOneComm(str, 0, cs, crs[0]);
  for(int ic= 1; ic< file->numComms; ++ic) {
    crs[ic]= crs[ic- 1]+ ncps;
    len= getline(&str, &n, file->fp);
    ncps= prvFile_readOneComm(str, ic, cs+ ic, crs[ic]);
  }
  if(NULL!= str) {
    free(str);
    str= NULL;
  }
  return 0;
}

inline static int PrvFile_reloadRecords(const ParaverFile *const file)
{
  int ret= fseeko(file->fp, file->recsPos, SEEK_SET);
  if(-1== ret) {
    fprintf(stderr, "%s: Cannot reload records!\n", __func__);
  }
  return ret;
}

inline static size_t prvFile_lastNewlinePos(const char *const buf,
                                            const size_t buflen, const size_t len)
{
  size_t ret= ULLONG_MAX;
  if('\0'== buf[0]|| 0== buflen) {
    char tmp[11]= { '\0' }; strncpy(tmp, buf, 10);
    fprintf(stderr, "%s: Empty buffer, returning ULLONG_MAX (buffer= \"%s\", len= %lu)\n", __func__, tmp, buflen);
    return ret;
  }
  size_t i= 0== len? buflen- 1: len- 1;
  for(; i> 0; --i) {
    if('\n'== buf[i]) {
      break;
    }
  }
  if('\n'== buf[i]) {
    ret= i;
  }
  return ret;
}
inline static void prvFile_processBuffer(char *const buf,
                                         void (*process)(char *const))
{
  char *ptr= strtok(buf, "\n");
  while(NULL!= ptr) {
    process(ptr);
    ptr= strtok(NULL, "\n");
  }
}
inline static double prvFile_timer_s()
{
  struct timespec ts;
  if(0!= clock_gettime(CLOCK_MONOTONIC, &ts)) {
    fprintf(stderr, "%s: Error obtaining clock-value\n", __func__);
    return -1.0;
  }
  return ts.tv_sec+ (ts.tv_nsec* 1.0e-9);
}

inline static double PrvFile_process(ParaverFile *const file,
                                     const int verbosity)
{
  if(NULL== file->lineProcessor) {
    return 0.0;
  }
  const size_t numBytes= (size_t) (file->size- file->recsPos);
  if(verbosity> 1) {
    if(verbosity> 2) {
      printf("Trace body size: %.1lf MB\n", ((double) numBytes)/ 1024.0/ 1024.0);
    }
    printf("Processed %02d%%...", 0); fflush(stdout);
  }
  const size_t buflen= 32* 1024* 1024;
  char *buf= malloc(sizeof(char)* (buflen+ 1)); buf[buflen]= '\0';
  size_t car= 0, rem= buflen- 1, numBytesRead= 0, numBytesProcessed= 0;

  FILE *fp= file->fp;
  double ioTime= -prvFile_timer_s();
  while(0!= (numBytesRead= fread(buf+ car, 1, rem, fp)+ car)) {
    ioTime+= prvFile_timer_s();
    size_t len= prvFile_lastNewlinePos(buf, buflen, numBytesRead);
    buf[len]= '\0';
    numBytesProcessed+= len+ 1;
    prvFile_processBuffer(buf, file->lineProcessor);
    car= numBytesRead- len- 1;
    rem= numBytesRead- car;
    memmove(buf, buf+ len+ 1, car);
    if(verbosity> 1) {
      printf("\rProcessed %02d%%...", (int) (numBytesProcessed* 100/ numBytes));
      fflush(stdout);
    }
    ioTime-= prvFile_timer_s();
  }
  ioTime+= prvFile_timer_s();
  if(verbosity> 1) {
    printf("\n"); fflush(stdout);
  }
  if(NULL!= buf) {
    free(buf);
    buf= NULL;
  }
  return ioTime;
}

/** @brief MPI function name table, indexed by Paraver event id (0-216). */
#define NUM_MPI_FUNCS 216
static const char *PrvMPINames[NUM_MPI_FUNCS]= {
  /* 0-8 */
  "Useful", "Send", "Recv", "Isend", "Irecv", "Wait", "Waitall", "Bcast", "Barrier",
  /* 9-15 */
  "Reduce", "Allreduce", "Alltoall", "Alltoallv", "Gather", "Gatherv", "Scatter",
  /* 16-21 */
  "Scatterv", "Allgather", "Allgatherv", "Comm_rank", "Comm_size", "Comm_create",
  /* 22-26 */
  "Comm_dup", "Comm_split", "Comm_group", "Comm_free", "Comm_remote_group",
  /* 27-31 */
  "Comm_remote_size", "Comm_test_inter", "Comm_compare", "Scan", "Init",
  /* 32-39 */
  "Finalize", "Bsend", "Ssend", "Rsend", "Ibsend", "Issend", "Irsend", "Test",
  /* 40-44 */
  "Cancel", "Sendrecv", "Sendrecv_replace", "Cart_create", "Cart_shift",
  /* 45-50 */
  "Cart_coords", "Cart_get", "Cart_map", "Cart_rank", "Cart_sub", "Cartdim_get",
  /* 51-55 */
  "Dims_create", "Graph_get", "Graph_map", "Graph_create", "Graph_neighbors",
  /* 56-60 */
  "Graphdims_get", "Graph_neighbors_count", "Topo_test", "Waitany", "Waitsome",
  /* 61-67 */
  "Probe", "Iprobe", "Win_create", "Win_free", "Put", "Get", "Accumulate",
  /* 68-73 */
  "Win_fence", "Win_start", "Win_complete", "Win_post", "Win_wait", "Win_test",
  /* 74-79 */
  "Win_lock", "Win_unlock", "Pack", "Unpack", "Op_create", "Op_free",
  /* 80-84 */
  "Reduce_scatter", "Attr_delete", "Attr_get", "Attr_put", "Group_difference",
  /* 85-89 */
  "Group_excl", "Group_free", "Group_incl", "Group_intersection", "Group_rank",
  /* 90-93 */
  "Group_range_excl", "Group_range_incl", "Group_size", "Group_translate_ranks",
  /* 94-97 */
  "Group_union", "Group_compare", "Intercomm_create", "Intercomm_merge",
  /* 98-102 */
  "Keyval_free", "Keyval_create", "Abort", "Error_class", "Errhandler_create",
  /* 103-106 */
  "Errhandler_free", "Errhandler_get", "Error_string", "Errhandler_set",
  /* 107-111 */
  "Get_processor_name", "Initialized", "Wtick", "Wtime", "Address",
  /* 112-116 */
  "Bsend_init", "Buffer_attach", "Buffer_detach", "Request_free", "Recv_init",
  /* 117-121 */
  "Send_init", "Get_count", "Get_elements", "Pack_size", "Rsend_init",
  /* 122-127 */
  "Ssend_init", "Start", "Startall", "Testall", "Testany", "Test_cancelled",
  /* 128-132 */
  "Testsome", "Type_commit", "Type_contiguous", "Type_extent", "Type_free",
  /* 133-137 */
  "Type_hindexed", "Type_hvector", "Type_indexed", "Type_lb", "Type_size",
  /* 138-143 */
  "Type_struct", "Type_ub", "Type_vector", "File_open", "File_close", "File_read",
  /* 144-147 */
  "File_read_all", "File_write", "File_write_all", "File_read_at",
  /* 148-151 */
  "File_read_at_all", "File_write_at", "File_write_at_all", "Comm_spawn",
  /* 152-155 */
  "Comm_spawn_multiple", "Request_get_status", "Ireduce", "Iallreduce",
  /* 156-161 */
  "Ibarrier", "Ibcast", "Ialltoall", "Ialltoallv", "Iallgather", "Iallgatherv",
  /* 162-167 */
  "Igather", "Igatherv", "Iscatter", "Iscatterv", "Ireducescat", "Iscan",
  /* 168-171 */
  "Reduce_scatter_block", "Ireduce_scatter_block", "Alltoallw", "Ialltoallw",
  /* 172-174 */
  "Get_accumulate", "Dist_graph_create", "Neighbor_allgather",
  /* 175-177 */
  "Ineighbor_allgather", "Neighbor_allgatherv", "Ineighbor_allgatherv",
  /* 178-180 */
  "Neighbor_alltoall", "Ineighbor_alltoall", "Neighbor_alltoallv",
  /* 181-183 */
  "Ineighbor_alltoallv", "Neighbor_alltoallw", "Ineighboralltoallw",
  /* 184-187 */
  "Fetch_and_op", "Compare_and_swap", "Win_flush", "Win_flush_all",
  /* 188-192 */
  "Win_flush_local", "Win_flush_local_all", "Mprobe", "Improbe", "Mrecv",
  /* 193-196 */
  "Imrecv", "Comm_split_type", "File_write_all_begin", "File_write_all_end",
  /* 197-199 */
  "File_read_all_begin", "File_read_all_end", "File_write_at_all_begin",
  /* 200-202 */
  "File_write_at_all_end", "File_read_at_all_begin", "File_read_at_alll_end",
  /* 203-205 */
  "File_read_ordered", "File_read_ordered_begin", "File_read_ordered_end",
  /* 206-208 */
  "File_read_shared", "File_write_ordered", "File_write_ordered_begin",
  /* 209-211 */
  "File_write_ordered_end", "File_write_shared", "Comm_dup_with_info",
  /* 212-215 */
  "Dist_graph_create_adjacent", "Comm_create_group", "Exscan", "Iexscan",
};
inline static const char *PrvFile_MPIName(const int eventId)
{
  if(eventId> -1&& eventId< NUM_MPI_FUNCS) {
    return PrvMPINames[eventId];
  }
  return NULL;
}
#undef NUM_MPI_FUNCS

/* -------------------------------------------------------------------------- */
/* END of private APIs of this library                                        */
/* -------------------------------------------------------------------------- */

#endif  /* CLOCKTALK_PARAVER_PARAVER_FILE_H__ */
