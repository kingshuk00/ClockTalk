# ClockTalk

- Trace replay for critical path from Paraver file.
- Usage:
  ```bash
  Usage: clocktalk [-PXT?V] [-m window,event] [-E[1]] [-R[1]] [--eager-limit=32k]
              [--ignore-events=traceability,flush,overhead]
              [--monitors=window,event] [--wmon-len=1.0e9] [--wmon-sma=1]
              [--emon-nevts=1] [--emon-rank=0] [--show-errors[=1]]
              [--show-reviews[=1]] [--pretty-output] [--export-profile]
              [--show-timings] [--help] [--usage] [--version]
              <paraver-file-name>
  ```
- ClockTalk works correctly only on MPI-traces.
- Calculates the length of the ideal-network constrined critical-path (aka "_ideal run time_") from Paraver traces.
  - assuming a network with _infinite bandwidth_ and _zero latency_
  - MPI communications starts <u>instantly</u> as soon as the communicating partners arrive - _zero latency_
  - When all partners have arrived, <u>communication finishes instantly irrespective of the message sizes</u> - _infinite bandwidth_
- The underlying methodology exchanges the Lamport clocks extended for MPI from sender to receiver (of messages) by asserting the causal relationship that the receive cannot finish (happen) before send, but the receive finishes as soon as the send is posted.
  - waits for the send to be posted (by receving the sender's Lamport clock).
  - starts a receive as soon as the process is in receive and the send is posted.
  - senders can choose between the following protocols:
    - `eager`: sending smaller messages finish right away without waiting for the receive relying on network hardware memory.
    - `rendezvous`: sends wait for matching recieve to be posted and then finish.
    - The default eager limit is `32 kB`.
- The model-factors monitoring framework creates a timeline of the model-factors for
  - fixed-windows including all ranks; or
  - event-driven window for one rank.
  - Monitors are an alternate to the timeline visualisation

## Obtain the code:
```bash
git clone git@github.com:kingshuk00/ClockTalk.git
git checkout main
```
## Compilation
- After checking out the `main` branch:
  ```bash
  cd ClockTalk
  scons -j8
  ```
- By default, the JSON compilation database for Clang tools is disabled. To enable,
  ```bash
  scons --compdb=1
  ```
- Default optimisation option is `-O3`. To compile with the debug symbols (`-g3 -O0`):
  ```bash
  scons --buildtype=debug
  ```

- Runtime options:
  - General usage:
    ```
    clocktalk --usage
    clocktalk --help
    clocktalk --version
    ```
  - Set diagnostic review level into `stdout`:
    ```bash
    clocktalk -R<n>
    clocktalk --show-reviews<=n>
    ```
  - Non-fatal logical errors encountered during processing of traces are reported in `stdout`. Turn them off by:
    ```
    clocktalk -E0
    clocktalk --show-errors=0
    ```
  - The default format of the end-output can be made pretty by:
    ```bash
    clocktalk -P
    clocktalk --pretty-output
    ```
  - After just reading the trace without replaying, quick profiling can be directed to a separate file.
    ```bash
    clocktalk -X
    clocktalk --export-profile
    ```
  - I/O progress and timings of steps can be displayed on `stdout`:
    ```bash
    clocktalk -T
    clocktalk --show-timings
    ```
  - Special trace events (trace-initialisation, flush, and trace-disable) are detected and discarded from calculations. These events can be ignored which results in such stretches treated as useful.
    ```bash
    clocktalk --ignore-events=[traceability,flush,overhead]
    ```
  - Default eager limit is 32kB. This can be changed:
    ```bash
    clocktalk --eager-limit=256k
    ```
    Default unit is `k`. Other valid units are `B`, `M`, `G`.
  - Use fixed-sized windowed monitoring:
    ```bash
      clocktalk -m window
    ```
    - Default window length is 1.0e9 ns. Change it with:
      ```bash
      clocktalk -m window --wmon-len 2.0e9
      ```
  - Use event-driven windowed monitoring:
    ```bash
      clocktalk -m event
    ```
    - Default rank is rank-0. Change it with:
      ```bash
      clocktalk -m event --emon-rank 1
      ```
      Be careful to use a valid rank.
    - Default number of events is 1. Change it with:
      ```bash
      clocktalk -m event --emon-nevts 8
      ```
  - The created file `<prv-filename.[w/e].dat` can be plotted using `gnuplot` with the script provided in the `utils` directory.
    ```gp
    gnuplot --persist -e "fname='<prv-filename>.[w/e]m.dat'" utils/mon.gp
    ```
    - Alternatively, Python can be used for creating the same plot (requires: `numpy`, `pandas`, `matplotlib`, `seaborn`)
      ```python
      python utils/mon.py <prv-filename>.[w/e]m.dat
      ```
