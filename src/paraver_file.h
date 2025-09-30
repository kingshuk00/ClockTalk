/*
 * Copyright (c) 2023-2025 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef PARAVER_FILE_READER_H__
#define PARAVER_FILE_READER_H__

#define _LARGEFILE_SOURCE

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<time.h>
#include<limits.h>
#include<assert.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>

typedef struct ParaverFile_struct__ {
  FILE *fp;
  off_t size;
  off_t communicatorsAt;
  off_t recordsAt;

  long long runTime;
  char timeUnit[4];
  int numNodes;
  int numApps;
  int numProcs;
  int numComms;
  long numAllCommsSizes;

  void (*lineProcessor)(char *const);
} ParaverFile;

inline static ParaverFile *ParaverFileOpen(const char *const filename);
inline static void ParaverFileClose(ParaverFile *const paraverFile);
inline static long long ParaverFileGetRuntime(const ParaverFile *const
                                              paraverFile) { return paraverFile->runTime; }
inline static const char *ParaverFileGetTimeUnit(const ParaverFile *const
                                                 paraverFile) { return paraverFile->timeUnit; }
inline static int ParaverFileGetNumNodes(const ParaverFile *const paraverFile) { return paraverFile->numNodes; }
inline static int ParaverFileGetNumApps(const ParaverFile *const paraverFile) { return paraverFile->numApps; }
inline static int ParaverFileGetNumProcs(const ParaverFile *const paraverFile) { return paraverFile->numProcs; };
inline static int ParaverFileGetNumComms(const ParaverFile *const paraverFile) { return paraverFile->numComms; }
inline static int ParaverFileGetAllCommsSizes(const ParaverFile *const
                                              paraverFile) { return paraverFile->numAllCommsSizes; }
inline static int ParaverFileReadComms(const ParaverFile *const paraverFile,
                                       int *const commsSizes,
                                       int **const commsRanks);

inline static int ParaverFileReloadRecords(const ParaverFile *const
                                           paraverFile);
inline static void ParaverFileSetLineProcessor(ParaverFile *const paraverFile,
                                               void (*processorFunc)(char *const))
{ paraverFile->lineProcessor= processorFunc; }
inline static double ParaverFileProcess(ParaverFile *const paraverFile,
                                         const bool silently);
inline static char *ParaverRecordNextNum(char *const p)
{
  char *x= strchr(p, ':');
  assert(NULL!= x);
  return x+ 1;
}
inline static char *ParaverRecordNextNumNth(char *const p, const int n)
{
  char *ptr= p;
  for(int i= 0; i< n; ++i) {
    ptr= ParaverRecordNextNum(ptr);
  }
  return ptr;
}

inline static const char *ParaverFileGetMPIName(const int eventId);

/* private */
inline static ParaverFile *ParaverFileOpen(const char *const filename)
{
  FILE *fp= NULL;
  char *headerStr= NULL;
  ParaverFile *file= (ParaverFile *) malloc(sizeof(ParaverFile));
  if(NULL== file) {
    printf("%s: Error allocating memroy\n", __func__);
    goto bad;
  }
  memset(file, 0, sizeof(ParaverFile));

  fp= fopen(filename, "r");
  if(NULL== fp) {
    printf("%s: Error opening file-\"%s\"\n", __func__, filename);
    goto bad;
  }

  {
    struct stat fprvStat;
    if(fstat(fileno(fp), &fprvStat) < 0) {
      printf("%s: Error reading file statistics\n", __func__);
      goto bad;
    }
    file->size= fprvStat.st_size;
  }

  if(-1== fseeko(fp, 0, SEEK_SET)) {
    printf("%s: Error locating Paraver header section.\n", __func__);
    goto bad;
  }

  size_t headerLen= 0;
  {
    ssize_t linelen= getline(&headerStr, &headerLen, fp);
    if(-1== linelen) {
      if(NULL!= headerStr) {
        free(headerStr);
        headerStr= NULL;
      }
      goto bad;
    } else {
      headerLen= (size_t) linelen;
    }
  }
  file->communicatorsAt= ftello(fp);

  char *ptr= headerStr;
  if(0!= strncmp("#Paraver (", ptr, 10)) {
    printf("%s: Unknown trace format - invalid header.\n", __func__);
  }
  ptr= strchr(ptr, ')')+ 2;
  file->runTime= atoll(ptr);
  if(0== strncmp("_ns", ptr- 3, 3)) {
    strcpy(file->timeUnit, "ns");
  } else {
    strcpy(file->timeUnit, "us");
  }

  ptr= strchr(ptr, ':')+ 1;
  file->numNodes= atoi(ptr);

  ptr= strchr(ptr, ':')+ 1;
  file->numApps= atoi(ptr);
  if(file->numApps> 1) {
    printf("Handling trac-files with more than 1 applications is not yet supported\n");
    printf("Results will be inaccurate.\n");
  }

  file->numProcs= 0;
  file->numComms= 0;
  ptr= strchr(ptr, ':')+ 1;
  {
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
  }

  /* printf("%d communicators found!\n", file->numComms); */
  for(int i= 0; i< file->numComms; ++i) {
    ssize_t len= getline(&headerStr, &headerLen, fp);
    if(-1== len) {
      printf("%s: Call to getline() failed.\n", __func__);
    }
    int cnp;
    /*                        c:app: id:np:p0:p1... */
    if(1!= sscanf(headerStr, "c:%*d:%*d:%d:", &cnp)) {
      printf("%s: Unknown trace format - cannot read communicators section.\n",
             __func__);
    }
    (file->numAllCommsSizes)+= cnp;
    /* printf("communicator-%d: %d members\n", id, cnp); */
  }
  file->recordsAt= ftello(fp);

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
  if(NULL!= headerStr) {
    free(headerStr);
    headerStr= NULL;
  }
  return file;
}
inline static void ParaverFileClose(ParaverFile *const file)
{
  if(NULL!= file) {
    if(NULL!= file->fp) {
      fclose(file->fp);
      file->fp= NULL;
    }
    free(file);
  }
}

inline static int paraverFileReadOneComm(char *const str, const int ix,
                                         int *const cs, int *const crs)
{
  int cix;
  /*                ignored-c:app: */
  if(2!= sscanf(str, "c:%*d:%d:%d:", &cix, cs)) {
    printf("%s: Unknown trace format - cannot read communicators section.\func",
           str);
  }
  --cix;
  if(cix!= ix) {
    printf("Problem, communicator-index erratic\n");
  }
  char *ptr= ParaverRecordNextNumNth(str, 4);
  for(int i= 0; i< *cs; ++i) {
    crs[i]= atoi(ptr)- 1;
    ptr= ParaverRecordNextNum(ptr);
  }
  return *cs;
}
inline static int ParaverFileReadComms(const ParaverFile *const file,
                                       int *const cs, int **const crs)
{
  if(file->numComms< 1) {
    return 0;
  }

  int ret= fseeko(file->fp, file->communicatorsAt, SEEK_SET);
  if(-1== ret) {
    printf("%s: cannot reload communicators!\n", __func__);
    return -1;
  }

  char *str; size_t n= 0;
  ssize_t len= getline(&str, &n, file->fp);
  if(-1== len) {
    printf("%s: Call to getline() failed.\n", __func__);
  }
  int ncps= paraverFileReadOneComm(str, 0, cs, crs[0]);
  for(int ic= 1; ic< file->numComms; ++ic) {
    crs[ic]= crs[ic- 1]+ ncps;
    len= getline(&str, &n, file->fp);
    ncps= paraverFileReadOneComm(str, ic, cs+ ic, crs[ic]);
  }
  if(NULL!= str) {
    free(str);
    str= NULL;
  }
  return 0;
}

inline static int ParaverFileReloadRecords(const ParaverFile *const file)
{
  int ret= fseeko(file->fp, file->recordsAt, SEEK_SET);
  if(-1== ret) {
    printf("%s: cannot reload records!\n", __func__);
  }
  return ret;
}

inline static size_t paraverFileGetLastNewlinePos(const char *const buf,
                                                  const size_t buflen, const size_t len)
{
  size_t ret= ULLONG_MAX;
  if('\0'== buf[0]|| 0== buflen) {
    char tmp[11]= { '\0' }; strncpy(tmp, buf, 10);
    printf("returning ULLONG_MAX (buffer= \"%s\", buflen= %lu\n", tmp, buflen);
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
inline static void paraverFileProcessBuffer(char *const buf,
                                            void (*process)(char *const))
{
  char *ptr= strtok(buf, "\n");
  while(NULL!= ptr) {
    process(ptr);
    ptr= strtok(NULL, "\n");
  }
}
inline static double IOTimer_s()
{
  struct timespec ts;
  if(0!= clock_gettime(CLOCK_MONOTONIC, &ts)) {
    printf("Error obtaining clock-value\n");
    return -1.0;
  }
  return ts.tv_sec+ (ts.tv_nsec* 1.0e-9);
}
/* returns time in seconds spent in fread() */
inline static double ParaverFileProcess(ParaverFile *const file,
                                        const bool silently)
{
  if(NULL== file->lineProcessor) {
    return 0.0;
  }
  const size_t numBytes= (size_t) (file->size- file->recordsAt);
  if(!silently) {
    printf("Size after comms section: %.1lf MB\n",
           ((double) numBytes)/ 1024.0/ 1024.0);
    printf("Processed %02d%%...", 0); fflush(stdout);
  }
  const size_t buflen= 32* 1024* 1024;
  char *buf= malloc(sizeof(char)* (buflen+ 1)); buf[buflen]= '\0';
  size_t car= 0, rem= buflen- 1, numBytesRead= 0, numBytesProcessed= 0;

  FILE *fp= file->fp;
  double ioTime= -IOTimer_s();
  while(0!= (numBytesRead= fread(buf+ car, 1, rem, fp)+ car)) {
    size_t len= paraverFileGetLastNewlinePos(buf, buflen, numBytesRead);
    buf[len]= '\0';
    numBytesProcessed+= len+ 1;
    ioTime+= IOTimer_s();
    paraverFileProcessBuffer(buf, file->lineProcessor);
    car= numBytesRead- len- 1;
    rem= numBytesRead- car;
    memmove(buf, buf+ len+ 1, car);
    if(!silently) {
      printf("\rProcessed %02d%%...", (int) (numBytesProcessed* 100/ numBytes));
      fflush(stdout);
    }
    ioTime-= IOTimer_s();
  }
  ioTime+= IOTimer_s();
  if(!silently) {
    printf("\n"); fflush(stdout);
  }
  if(NULL!= buf) {
    free(buf);
    buf= NULL;
  }
  return ioTime;
}

#define NUM_MPI_FUNCS 194
static const char *ParaverMPINames[NUM_MPI_FUNCS]= {
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
  /* 193-193 */
  "Imrecv",
};
inline static const char *ParaverFileGetMPIName(const int ev)
{
  if(ev> -1&& ev< NUM_MPI_FUNCS) {
    return ParaverMPINames[ev];
  }
  return NULL;
}
#undef NUM_MPI_FUNCS

#endif  /* PARAVER_FILE_READER_H__ */
