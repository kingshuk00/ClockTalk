/*
 * Copyright (c) 2026      Kingshuk Haldar. All rights reserved.
 *
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart. All rights reserved.
 *
 * Authors: Kingshuk Haldar <haldar.kingshuk@gmail.com>
 *
 */


#include<stdint.h>

/**
 * @file paraver_file.h
 * @brief General-purpose header-only library to aid read Paraver trace files.
 *
 * A Paraver trace file has three sections in order:
 *
 *   1. **Header line**: a single line beginning with `#Paraver (` that encodes
 *      the trace duration, time unit, node count, CPUs per node, application
 *      count, process counts per application, and communicator counts.
 *
 *   2. **Communicator section**: one line description per communicator, each
 *      in the form
 *      `c:<app-id>:<comm-id>:<comm-size>:<comm-rank-0>:<comm-rank-1>:...`,
 *      listing the global MPI ranks that belong to each communicator.
 *
 *   3. **Records section**: main body of the trace. Every line begins with a
 *      record-type digit: `1: state`, `2: event`, `3: point-to-point
 *      message`.
 *
 * `ParaverFile_open()` function parses the header and communicator section
 * entirely and stores the hardware topology, application layout, and
 * communicator information within the `ParaverFile` struct.
 *
 * Afterwards, the consumer of this library calls `ParaverFile_process()` to
 * process through the records section, one record-line at a time via a
 * callback set with the `ParaverFile_setLineProcessor()` function.
 * A multi-pass workflow is supported by calling `ParaverFile_reloadRecords()`
 * between passes to seek back to the beginning of the records section.
 *
 * The records section is read in `32 MB` chunks using `fread()` for throughput.
 * Partial lines at chunk boundaries are carried over to the next chunk.
 * The time spent in `fread()` is measured separately from processing
 * time and returned by the `ParaverFile_process()` function to the consumer.
 *
 * @note
 * - Requires `_LARGEFILE_SOURCE` for `fseeko()`/`ftello()` on 32-bit systems
 *   so that file offsets are 64-bit and traces larger than 2 GB are handled
 *   correctly. On 64-bit systems this has no effect.
 *
 * - Multi-application traces are not supported. `ParaverFile_open()` issues
 *   a fatal error in such scenarios.
 *
 * - Designed to be used standalone without any external dependencies beyond
 *   the C and POSIX standards.
 *
 * - Consumers must define a label `bad` in function for exiting / handling
 *   of errors when `prv_err()` function is called from a public API.
 */

#ifndef CLOCKTALK_PARAVER_PARAVER_FILE_H__
#define CLOCKTALK_PARAVER_PARAVER_FILE_H__

#define _LARGEFILE_SOURCE

/******************************************************************************/
/* Public APIs of this header-only library                                    */
/******************************************************************************/

typedef struct ParaverFile_struct__ ParaverFile;

/**
 * @brief Opens a Paraver trace file and parses its header and communicator
 *        sections.
 *
 * Allocates and returns a `ParaverFile` handle populated with complete metadata
 * from the header and communicator sections. The file pointer is positioned at
 * the start of the records section on success.
 *
 * Specifically:
 *   - Opens the file and records its total size via `fstat()`.
 *   - Reads and parses the header and communicators sections to extract
 *     runtime, time unit, node count, per-node CPU counts, application count,
 *     per-application task counts, thread/node assignments, and communicator
 *     counts.
 *   - Records file offsets for the communicator and records sections so they
 *     can be seeked to later.
 *
 * @param[in] fn Path to the `.prv` trace file.
 * @return Allocated and initialised `ParaverFile` handle, or `NULL` on an error
 *         (eg. file not found, unrecognised header format, allocation failure).
 *         All resources are cleaned up before returning `NULL`.
 */
inline static ParaverFile *ParaverFile_open(const char *const fn);

/**
 * @brief Forgets the metadata structs in `file`.
 *
 * Sets the application, hardware topology, communicator structs as `NULL`.
 * Consumers can call this after copying those addresses to avoid them getting
 * `free()`d when `ParaverFile_close()` is called.
 *
 * @param[in] file Path to the `.prv` trace file.
 */
inline static void ParaverFile_disownMetadata(ParaverFile *const file);

/**
 * @brief Closes the file handle and deep-frees the `ParaverFile` struct.
 *
 * Safe to call with `NULL`.
 *
 * @param[in] file Handle to close. Must not be used after this call.
 */
inline static void ParaverFile_close(ParaverFile *const file);

/**
 * @brief Returns the trace duration in the trace's native time unit.
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return Trace duration as parsed from the header.
 */
inline static int64_t ParaverFile_duration(const ParaverFile *const file);

/**
 * @brief Returns the time unit string declared in the trace header.
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return Time unit (e.g. "ns") as a pointer to an immutable NULL-terminated C
 *         string. If the header does not declare anything, "us" is returned.
 */
inline static const char *ParaverFile_timeUnit(const ParaverFile *const file);

/**
 * @brief Returns the number of hardware nodes recorded in the header.
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return Node count.
 */
inline static int32_t ParaverFile_numNodes(const ParaverFile *const file);

/**
 * @brief Returns the number of CPUs on a specific hardware node.
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @param[in] node 0-based node index. Must be in [0, `ParaverFile_numNodes()`).
 * @return CPU count for the given node as recorded in the header.
 */
inline static int32_t ParaverFile_numNodeCPUs(const ParaverFile *const file, const int32_t node);

/**
 * @brief Returns the total number of CPUs across all hardware nodes.
 *
 * Equivalent to summing `ParaverFile_numNodeCPUs()` over all nodes.
 *
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return Total CPU count across all nodes.
 */
inline static int32_t ParaverFile_numCPUs(const ParaverFile *const file);

/**
 * @brief Returns the number of applications recorded in the header.
 *
 * Multi-application traces are treated as a fatal error by `ParaverFile_open()`
 * and will never be seen here in normal usage.
 *
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return Application count, always 1 for supported traces.
 */
inline static int32_t ParaverFile_numApps(const ParaverFile *const file);

/**
 * @brief Returns the number of tasks (MPI ranks) in a specific application.
 *
 * In Paraver terminology, a task corresponds to one MPI process.
 *
 * @note Threads within a task are tracked separately in the header but not used
 * by ClockTalk, which operates at the task (rank) granularity.
 *
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @param[in] app  0-based application index. Must be in
 *                 [0, `ParaverFile_numApps()`).
 * @return Number of tasks in the given application.
 */
inline static int32_t ParaverFile_numAppTasks(const ParaverFile *const file, const int32_t app);

/**
 * @brief Returns the total number of tasks (MPI ranks) across all applications.
 *
 * Equivalent to summing `ParaverFile_numAppTasks()` over all applications.
 * This is the value to use when allocating per-process arrays.
 *
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return Total task count across all applications.
 */
inline static int32_t ParaverFile_numTasks(const ParaverFile *const file);
/**
 * @brief Returns the number of communicators declared in the trace.
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return Communicator count.
 */
inline static int32_t ParaverFile_numComms(const ParaverFile *const file);

/**
 * @brief Returns the number of member ranks in a specific communicator.
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @param[in] comm 0-based communicator index. Must be in
 *                 [0, `ParaverFile_numComms()`).
 * @return Number of MPI ranks in the communicator.
 */
inline static int32_t ParaverFile_commSize(const ParaverFile *const file, const int32_t comm);

/**
 * @brief Returns the global MPI rank of the i-th member of a communicator.
 *
 * Ranks are 0-based global task indices as used throughout the trace data.
 *
 * @param[in] file  Opened file handle returned by `ParaverFile_open()`.
 * @param[in] comm  0-based communicator index.
 *                  Must be in [0, `ParaverFile_numComms()`).
 * @param[in] irank 0-based member index within the communicator.
 *                  Must be in [0, `ParaverFile_commSize(file, comm)`).
 * @return 0-based global MPI rank of the given member.
 */
inline static int32_t ParaverFile_commRank(const ParaverFile *const file, const int32_t comm, const int32_t irank);

/**
 * @brief Seeks the file pointer back to the start of the records section.
 *
 * This function can be used in-between two passes at the consumer (e.g. a count
 * pass and a read pass) to restart record processing without reopening the
 * file.
 *
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @return 0 on success, -1 on seek failure.
 */
inline static int ParaverFile_reloadRecords(const ParaverFile *const file);

/**
 * @brief Sets the line callback invoked for each record line during processing.
 *
 * The callback receives a null-terminated string containing one record line
 * with the trailing newline replaced by '\0'. It is called once per line
 * during `ParaverFile_process()`.
 *
 * Must be called before `ParaverFile_process()` to set the callback.
 * Call it again with a different function between passes.
 *
 * @param[in] file Opened file handle returned by `ParaverFile_open()`.
 * @param[in] processorFunc Callback function that processes one mutable record
 *                          line at a time. May modify the string in place
 *                          (e.g. via `strtok()`) but must not free it.
 */
inline static void ParaverFile_setLineProcessor(ParaverFile *const file,
                                                void (*processorFunc)(char *const));

/**
 * @brief Reads and processes the entire records section.
 *
 * Reads the file in 32 MB chunks using `fread()` for throughput. Partial lines
 * at chunk boundaries are carried forward with `memmove()`. The line processor
 * set by `ParaverFile_setLineProcessor()` is called once for each complete line.
 *
 * The time spent inside `fread()` is measured separately from processing time
 * using `CLOCK_MONOTONIC` and returned to the caller for I/O performance
 * reporting.
 *
 * @param[in] file Opened file handle returned by `ParaverFile_open()`. Must have a
                   line processor set via `ParaverFile_setLineProcessor()`.
 * @param[in] verbosity Controls progress output to `stdout`. `0`: silent;
 *                      non-zero: progress percentage to `stdout`.
 * @return Time spent in `fread()` in seconds. Returns 0.0 if no line processor
 *         has been set.
 */
inline static double ParaverFile_process(ParaverFile *const file,
                                         const int verbosity);

/**
 * @brief Advances a pointer past the next occurrence of delimiter `c`.
 *
 * Paraver record lines use `:` as the primary field delimiter, but the header
 * also uses `(`, `)`, and `,` as delimiters. This generalised version handles
 * any delimiter character.
 *
 * Returns a pointer to the character immediately after the first occurrence of
 * `c` in `p`, i.e. to the start of the next field value. Returns `NULL` if
 * `c` is not found, allowing the caller to detect end-of-line cleanly without
 * any potential assertion failure.
 *
 * @param[in] p Pointer into a record or header line.
 * @param[in] c Delimiter character to search for.
 * @return Pointer to the first character after `c`, or `NULL` if not found.
 */
inline static const char *ParaverFile_recAfter(const char *const p, const char c);

/**
 * @brief Advances a pointer past `n` consecutive occurrences of delimiter `c`.
 *
 * Calls `ParaverFile_recAfter()` repeatedly `n` times. Returns `NULL` as soon
 * as any intermediate call returns `NULL`, so the caller can detect a
 * truncated or malformed line without further checks.
 *
 * @param[in] p Pointer into a record or header line.
 * @param[in] c Delimiter character to skip past.
 * @param[in] n Number of delimiter occurrences to skip.
 * @return Pointer to the first character after the nth occurrence of `c`,
 *         or `NULL` if fewer than `n` occurrences exist.
 */
inline static const char *ParaverFile_recAfterNth(const char *const p, const char c, const int n);

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
inline static const char *ParaverFile_MPIName(const int eventId);

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
 * Populated entirely by `ParaverFile_open()`, which parses the header and the
 * communicator sections.
 *
 * @note All fields should be considered read-only. Exceptions are the following
  - `ParaverFile_disownMetadata()` that releases the ownership; and
  - `ParaverFile_setLineProcessor()` that is mutated during usage.
*/
typedef struct ParaverFile_struct__ {
  FILE *fp;
  off_t size;

  off_t commsPos;
  off_t recsPos;

  int verbosity;

  struct {
    /** @brief Trace duration from the trace file */
    int64_t duration;

    /** @brief NULL-terminated string for the duration unit */
    char unit[4];
  } runtime;

  struct {
    /** @brief Number of CPUs on each node. */
    int32_t numCPUs;
  } *nodes;
  int32_t numNodes;

  struct {
    /** @brief Array of `numTasks` task descriptors. */
    struct {
      /** @brief Number of threads in each MPI ranks. */
      int32_t numThreads;

      /** @brief Node-id where this task belongs. */
      int32_t nodeId;
    } *tasks;

    /** @brief Number of tasks (MPI ranks) in each application. */
    int32_t numTasks;
  } *apps;
  int32_t numApps;

  struct {
    /** @brief malloc-ed array of global MPI ranks in this comm. */
    int32_t *ranks;

    /** @brief Number of ranks in this comm. */
    int32_t size;

    /** @brief Index of application to which this communicator belongs. */
    int32_t app;

    /** @brief Communicator index as read from the file. */
    int32_t ix;
  } *comms;
  int32_t numComms;

  void (*lineProcessor)(char *const);
} ParaverFile;


/**
 * @var ParaverFile::fp
 * @brief Opened Paraver file handle.
 * @details Positioned at the start of the records section after
 *          `ParaverFile_open()` returns, then repositioned by
 *          `ParaverFile_reloadRecords()`.
 *
 * @var ParaverFile::size
 * @brief Total file size in bytes, obtained via `fstat()`.
 * @details Used by `ParaverFile_process()` to compute and display progress.
 *
 * @var ParaverFile::commsPos
 * @brief Trace file offset of position for the communicator section.
 *
 * @var ParaverFile::recsPos
 * @brief Trace file offset of position for the record section.
 * @details `ParaverFile_reloadRecords()` resets current position of opened file
 *          here for multi-pass processing.
 *
 * @var ParaverFile::verbosity
 * @brief Verbosity level forwarded to `ParaverFile_process()`.
 *
 * @var ParaverFile::runTime
 * @brief Trace duration and its time unit.
 *
 * @var ParaverFile::nodes
 * @brief Array of `numNodes` node descriptors.
 *
 * @var ParaverFile::numNodes
 * @brief Number of hardware nodes recorded in the header
 *
 * @var ParverFile::apps
 * @brief Array of `numApps` application descriptors.
 *
 * @var ParaverFile::numApps
 * @brief Number of applications from the header.
 * @details Always 1 for supported traces; `ParaverFile_open()` treats
 *          numApps > 1 as a fatal error.
 *
 * @var ParaverFile::numComms
 * @brief Number of communicators declared in the trace.
 *
 * @var ParaverFile::lineProcessor
 * @brief Callback set by `ParaverFile_setLineProcessor()`.
 * @details Called once per record line during `ParaverFileProcess()`.
 */

inline static int64_t ParaverFile_duration(const ParaverFile *const file)
{
  return file->runtime.duration;
}
inline static const char *ParaverFile_timeUnit(const ParaverFile *const file)
{
  return file->runtime.unit;
}
inline static int32_t ParaverFile_numNodes(const ParaverFile *const file)
{
  return file->numNodes;
}
inline static int32_t ParaverFile_numNodeCPUs(const ParaverFile *const file, const int32_t node)
{
  return file->nodes[node].numCPUs;
}
inline static int32_t ParaverFile_numCPUs(const ParaverFile *const file)
{
  int32_t ncpus= 0;
  for(int32_t n= 0; n< ParaverFile_numNodes(file); ++n) {
    ncpus+= ParaverFile_numNodeCPUs(file, n);
  }
  return ncpus;
}
inline static int32_t ParaverFile_numApps(const ParaverFile *const file)
{
  return file->numApps;
}
inline static int32_t ParaverFile_numAppTasks(const ParaverFile *const file, const int32_t app)
{
  return file->apps[app].numTasks; /* task is MPI rank */
}
inline static int32_t ParaverFile_numTasks(const ParaverFile *const file)
{
  int32_t ntasks= 0;
  for(int32_t a= 0; a< ParaverFile_numApps(file); ++a) {
    ntasks+= ParaverFile_numAppTasks(file, a);
  }
  return ntasks;                /* task is MPI rank */
}
inline static int32_t ParaverFile_numComms(const ParaverFile *const file)
{
  return file->numComms;
}
inline static int32_t ParaverFile_commSize(const ParaverFile *const file, const int32_t comm)
{
  return file->comms[comm].size;
}
inline static int32_t ParaverFile_commRank(const ParaverFile *const file, const int32_t comm, const int32_t irank)
{
  return file->comms[comm].ranks[irank];
}

inline static void ParaverFile_setLineProcessor(ParaverFile *const file,
                                                void (*processorFunc)(char *const))
{
  file->lineProcessor= processorFunc;
}

inline static const char *ParaverFile_recAfter(const char *const p, const char c)
{
  const char *x= strchr(p, c);
  return NULL== x? x: x+ 1;
}

inline static const char *ParaverFile_recAfterNth(const char *const p, const char c, const int n)
{
  const char *ptr= p;
  for(int i= 0; i< n; ++i) {
    if(NULL== ptr) {
      break;
    }
    ptr= ParaverFile_recAfter(ptr, c);
  }
  return ptr;
}

#define paraver_error(...) do {                                         \
    fprintf(stderr, "*** "__VA_ARGS__);                                 \
    fprintf(stderr, "\n");                                              \
    fflush(stderr);                                                     \
    goto bad;                                                           \
  } while(0)

inline static off_t paraverfile_getSize(FILE *fp)
{
  struct stat fpStat;
  if(fstat(fileno(fp), &fpStat)< 0) {
    paraver_error("paraver-file-stat: failed fpstat() call");
  }

  if(-1== fseeko(fp, 0, SEEK_SET)) {
    paraver_error("paraver-file-seek: failed fseeko() call");
  }

  goto bye;

 bad:
  return -1;

 bye:
  return fpStat.st_size;
}

inline static ParaverFile *paraverfile_freeAll(ParaverFile *const file)
{
  if(NULL== file) {
    return NULL;
  }

  if(NULL!= file->comms) {
    for(int c= 0; c< file->numComms; ++c) {
      free(file->comms[c].ranks);
    }
    free(file->comms);
  }

  if(NULL!= file->apps) {
    for(int a= 0; a< file->numApps; ++a) {
      free(file->apps[a].tasks);
    }
    free(file->apps);
  }

  free(file->nodes);

  free(file);

  return NULL;
}

inline static ParaverFile *ParaverFile_open(const char *const fn)
{
  char *line= NULL;
  FILE *fp= NULL;
  ParaverFile *file= (ParaverFile *) malloc(sizeof(ParaverFile));
  if(NULL== file) {
    paraver_error("paraver-file-alloc: failed to allocate memroy for Paraver file struct");
  }
  memset(file, 0, sizeof(ParaverFile));

  fp= fopen(fn, "r");
  if(NULL== fp) {
    paraver_error("paraver-file-open: failed to open file-\"%s\"", fn);
  }

  if((file->size= paraverfile_getSize(fp))< 0) {
    paraver_error("paraver-file-size: failed to obtain size of opend file-\"%s\"", fn);
  }

  size_t len= 0;
  ssize_t slen= getline(&line, &len, fp);
  if(-1== slen) {
    paraver_error("paraver-file-header: failed getline() call on Paraver file-\"%s\"", fn);
  }
  file->commsPos= ftello(fp);

  const char *ptr= line;
  if(0!= strncmp("#Paraver (", ptr, 10)) {
    paraver_error("paraver-file-format: invalid beginning of Paraver file header");
  }

  ptr= ParaverFile_recAfter(ptr, ')');
  if(NULL== ptr) {
    paraver_error("paraver-file-format: invalid date format of Paraver file header");
  }

  ptr= ParaverFile_recAfter(ptr, ':');
  if(NULL== ptr) {
    paraver_error("paraver-file-format: invalid beginning of Paraver file duration");
  }
  file->runtime.duration= (int64_t) atoll(ptr);

  ptr= ParaverFile_recAfter(ptr, ':');
  if(NULL== ptr) {
    paraver_error("paraver-file-format: invalid beginning of Paraver file node information");
  }
  if(0== strncmp("_ns", ptr- 4, 3)) {
    strcpy(file->runtime.unit, "ns");
  } else {
    strcpy(file->runtime.unit, "us");
  }

  file->numNodes= (int32_t) atoi(ptr);
  ptr= ParaverFile_recAfter(ptr, '(');
  if(NULL== ptr) {
    paraver_error("paraver-file-format: invalid Paraver file node information");
  }
  file->nodes= (typeof(file->nodes)) malloc(sizeof(*(file->nodes))* file->numNodes);
  if(NULL== file->nodes) {
    paraver_error("paraver-file-alloc: failed to allocate memory for Paraver file nodes.");
  }
  for(int n= 0; n< file->numNodes- 1; ++n) {
    if(NULL== ptr) {
      paraver_error("paraver-file-format: invalid Paraver file node information");
    }
    file->nodes[n].numCPUs= atoi(ptr);
    ptr= ParaverFile_recAfter(ptr, ',');
  }
  if(NULL== ptr) {
    paraver_error("paraver-file-format: invalid Paraver file node information");
  }
  file->nodes[file->numNodes- 1].numCPUs= atoi(ptr);
  ptr= ParaverFile_recAfter(ptr, ')');
  if(NULL== ptr) {
    paraver_error("paraver-file-format: invalid end of Paraver file node information");
  }

  ptr= ParaverFile_recAfter(ptr, ':');
  if(NULL== ptr) {
    paraver_error("paraver-file-format: invalid Paraver file application count");
  }
  file->numApps= atoi(ptr);
  if(file->numApps> 1) {
    paraver_error("paraver-file-format: multi-application Paraver traces are not "
            "supported.\n                      Results will be inaccurate.");
  }

  file->apps= (typeof(file->apps)) malloc(sizeof(*(file->apps))* file->numApps);
  if(NULL== file->apps) {
    paraver_error("paraver-file-alloc: failed to allocate memory for Paraver file applications");
  }
  for(uint16_t a= 0; a< file->numApps; ++a) {
    ptr= ParaverFile_recAfter(ptr, ':');
    if(NULL== ptr) {
      paraver_error("paraver-file-format: invalid Paraver file application task count");
    }

    file->apps[a].numTasks= atoi(ptr);
    file->apps[a].tasks= (typeof(file->apps[a].tasks)) malloc(sizeof(*(file->apps[a].tasks))* file->apps[a].numTasks);
    if(NULL== file->apps[a].tasks) {
      paraver_error("paraver-file-alloc: failed to allocate memory for Paraver file tasks");
    }
    ptr= ParaverFile_recAfter(ptr, '(');
    if(NULL== ptr) {
      paraver_error("paraver-file-format: invalid Paraver file application tasks");
    }
    for(int t= 0; t< file->apps[a].numTasks- 1; ++t) {
      file->apps[a].tasks[t].numThreads= atoi(ptr);
      ptr= ParaverFile_recAfter(ptr, ':');
      if(NULL== ptr) {
        paraver_error("paraver-file-format: invalid Paraver file application tasks");
      }
      file->apps[a].tasks[t].nodeId= atoi(ptr);
      ptr= ParaverFile_recAfter(ptr, ',');
      if(NULL== ptr) {
        paraver_error("paraver-file-format: invalid Paraver file application tasks");
      }
    }
    file->apps[a].tasks[file->apps[a].numTasks- 1].numThreads= atoi(ptr);
    ptr= ParaverFile_recAfter(ptr, ':');
    if(NULL== ptr) {
      paraver_error("paraver-file-format: invalid Paraver file application tasks");
    }
    file->apps[a].tasks[file->apps[a].numTasks- 1].nodeId= atoi(ptr);
    ptr= ParaverFile_recAfter(ptr, ')');
    if(NULL== ptr) {
      paraver_error("paraver-file-format: invalid Paraver file application tasks");
    }

    /* expect ',' between *and* at the end of apps */
    /* <extrae/merger/paraver/paraver_generator.c:971 Paraver_WriteHeader() */
    ptr= ParaverFile_recAfter(ptr, ',');
    if(NULL== ptr) {
      paraver_error("paraver-file-format: invalid Paraver file application tasks");
    }
  }

  file->numComms= atoi(ptr);
  file->comms= (typeof(file->comms)) malloc(sizeof(*(file->comms))* file->numComms);
  if(NULL== file->comms) {
    paraver_error("paraver-file-alloc: failed to allocate memory for Paraver file communicators");
  }

  for(int c= 0; c< file->numComms; ++c) {
    slen= getline(&line, &len, fp);
    if(-1== slen) {
      paraver_error("paraver-file-format: failed reading Paraver file communicators");
    }
    /* c:app:comm:comm-size:comm-ranks... */
    if('c'!= line[0]) {
      paraver_error("paraver-file-format: invalid Paraver file communicator format");
    }
    ptr= ParaverFile_recAfter(line, ':');

    file->comms[c].app= atoi(ptr)- 1;
    ptr= ParaverFile_recAfter(ptr, ':');

    file->comms[c].ix= atoi(ptr)- 1;
    if(file->comms[c].ix!= c) {
      paraver_error("paraver-file-format: erratic communicator index (%d!= exp:%d)",
              file->comms[c].ix, c);
    }
    ptr= ParaverFile_recAfter(ptr, ':');

    file->comms[c].size= atoi(ptr);
    file->comms[c].ranks= (int32_t *) malloc(sizeof(int32_t)* file->comms[c].size);
    memset(file->comms[c].ranks, 0, sizeof(int32_t)* file->comms[c].size);
    for(int r= 0; r< file->comms[c].size; ++r) {
      ptr= ParaverFile_recAfter(ptr, ':');
      file->comms[c].ranks[r]= atoi(ptr)- 1;
    }
  }
  file->recsPos= ftello(fp);

  file->fp= fp;
  goto bye;

bad:
  if(NULL!= fp) {
    fclose(fp);
    fp= NULL;
  }

  file= paraverfile_freeAll(file);

bye:
  if(NULL!= line) {
    free(line);
    line= NULL;
  }

  return file;
}

inline static void ParaverFile_disownMetadata(ParaverFile *const file)
{
  file->nodes= NULL;
  file->numNodes= 0;

  file->apps= NULL;
  file->numApps= 0;

  file->comms= NULL;
  file->numComms= 0;
}

inline static void ParaverFile_close(ParaverFile *const file)
{
  if(NULL!= file) {
    if(NULL!= file->fp) {
      fclose(file->fp);
      file->fp= NULL;
    }
    paraverfile_freeAll(file);
  }
}

inline static int ParaverFile_reloadRecords(const ParaverFile *const file)
{
  int ret= fseeko(file->fp, file->recsPos, SEEK_SET);
  if(-1== ret) {
    fprintf(stderr, "%s: Cannot reload records!\n", __func__);
  }
  return ret;
}

inline static size_t paraverFile_lastNewlinePos(const char *const buf,
                                                const size_t buflen,
                                                const size_t len)
{
  size_t ret= ULLONG_MAX;
  if('\0'== buf[0]|| 0== buflen) {
    char tmp[11]= { '\0' }; strncpy(tmp, buf, 10);
    fprintf(stderr, "%s: Empty buffer, returning ULLONG_MAX "
            "(buffer= \"%s\", len= %lu)\n", __func__, tmp, buflen);
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
inline static void paraverFile_processBuffer(char *const buf,
                                             void (*process)(char *const))
{
  char *ptr= strtok(buf, "\n");
  while(NULL!= ptr) {
    process(ptr);
    ptr= strtok(NULL, "\n");
  }
}
inline static double paraverFile_timer_s()
{
  struct timespec ts;
  if(0!= clock_gettime(CLOCK_MONOTONIC, &ts)) {
    fprintf(stderr, "%s: Error obtaining clock-value\n", __func__);
    return -1.0;
  }
  return ts.tv_sec+ (ts.tv_nsec* 1.0e-9);
}

inline static double ParaverFile_process(ParaverFile *const file,
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
  double ioTime= -paraverFile_timer_s();
  while(0!= (numBytesRead= fread(buf+ car, 1, rem, fp)+ car)) {
    ioTime+= paraverFile_timer_s();
    size_t len= paraverFile_lastNewlinePos(buf, buflen, numBytesRead);
    buf[len]= '\0';
    numBytesProcessed+= len+ 1;
    paraverFile_processBuffer(buf, file->lineProcessor);
    car= numBytesRead- len- 1;
    rem= numBytesRead- car;
    memmove(buf, buf+ len+ 1, car);
    if(verbosity> 1) {
      printf("\rProcessed %02d%%...", (int) (numBytesProcessed* 100/ numBytes));
      fflush(stdout);
    }
    ioTime-= paraverFile_timer_s();
  }
  ioTime+= paraverFile_timer_s();
  if(verbosity> 1) {
    printf("\n"); fflush(stdout);
  }
  if(NULL!= buf) {
    free(buf);
    buf= NULL;
  }
  return ioTime;
}

#define NUM_PARAVER_MPI_FUNCS 216
/** @brief MPI function name table, indexed by Paraver event id (0-216). */
static const char *ParaverMPINames[NUM_PARAVER_MPI_FUNCS]= {
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
inline static const char *ParaverFile_MPIName(const int eventId)
{
  if(eventId> -1&& eventId< NUM_PARAVER_MPI_FUNCS) {
    return ParaverMPINames[eventId];
  }
  return NULL;
}
#undef NUM_PARAVER_MPI_FUNCS
#undef paraver_error

/* -------------------------------------------------------------------------- */
/* END of private APIs of this library                                        */
/* -------------------------------------------------------------------------- */

#endif  /* CLOCKTALK_PARAVER_PARAVER_FILE_H__ */
