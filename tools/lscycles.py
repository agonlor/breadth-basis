#!/usr/bin/python3

from os import sep
import sys
import json
import string

class Solutions():
  def __init__(self, args):
    # Initialize the arguments
    self.filenames = [] # filenames of the solutions
    self.sheet = False  # show the results for the spreadsheet
    self.average = False
    self.algorithms = ["tb", "deylw_2018"]
    
    for arg in args:
      if arg.startswith("-"):
        if "sheet" in arg:
          self.sheet = True
        if "average" in arg:
          self.average = True
      else:
        self.filenames.append(arg)

    # print usage if necessary
    if len(self.filenames) == 0:
      print("""Usage: ./lscycles [options] files
    Options:
      -sheet    Produce output formatted for a spreadsheet
      -average  Compute the average of solutions computed with the same algorithm
    Examples:
      - List all solutions in the current directory:  ./lscycles *.json
      - Produce a spreadsheet column:                 ./lscycles -sheet *.json
      - Print average statistics of each algorithm:   ./lscycles -average *.json""")
      exit(1)


  def read_solutions(self):
    """
    Read the solutions and save some statistics
    """
    self.instance_output = {} # map instance -> list of solutions
    for fn in self.filenames:
      try:
        with open(fn) as sol_file:
          sol = json.load(sol_file)
        stats = {} # information about the current solution
        stats['fn'] = fn
        stats['object'] = sol['meta']['object']
        stats['cycles'] = [len(cycle) for cycle in sol["cycles"]]
        stats['length'] = sum(stats['cycles'])
        stats['algo'] = sol['meta']['algorithm']
        stats['args'] = sol['meta']['arguments']
        stats['time'] = sol['meta']['elapsed_time']
        # initialize the list of solutions for this instance, if necessary
        if stats['object'] not in self.instance_output:
          self.instance_output[stats['object']] = []  
        self.instance_output[stats['object']].append(stats)
      except:
        print(f"Error reading solution file {fn}")
  

  def make_line(self, sol, fnlen):
    if self.sheet:
      line = f"{sol['fn']};{sol['length']};{sol['algo']};"
      line += "_".join(f"{arg}" for arg in sol['args'])
    else:
      line = f"{sol['fn']:{fnlen}} {sol['length']:4} {sol['algo']:11} "
      line += "_".join(f"{arg:2}" for arg in sol['args'])
      line += f" {sol['time']:.2f}"
      line += f" {sol['cycles']}"
    return line

  def print(self):
    if self.average:
      self.print_average()
    else:
      self.print_all()


  def print_all(self):
    fnlen = max(len(sol['fn']) for sols in self.instance_output.values() for sol in sols) # maximum length of the filenames
    fnlen += 2 # extra padding
    outs = ""

    for sols in self.instance_output.values():
      for sol in sols:
        outs += self.make_line(sol, fnlen) + "\n"
    print(outs)
  

  def print_average(self):
    """ Print the average result of each algorothm on each instance """
    if self.sheet:
      print("filename;" + ";".join(f"{algo}" for algo in self.algorithms))

    fnlen = max(len(obj) for obj in self.instance_output) # maximum length of the filenames
    for obj in self.instance_output:
      sols = {}
      for sol in self.instance_output[obj]:
        if sol['algo'] not in sols:
          sols[sol['algo']] = []
        sols[sol['algo']].append(sol['length'])
      if self.sheet:
        print(f'{obj:{fnlen}}', end=';')
        for algo in self.algorithms:
          print(f"{sum(sols[algo])/len(sols[algo]):.2f}", end=';')
        print("")
      else:
        print(f'{obj:{fnlen}}:', end=' ')
        for algo in self.algorithms:
          print(f"{algo}={sum(sols[algo])/len(sols[algo]):.2f} ({len(sols[algo])})", end=' ')
        print("")


solutions = Solutions(sys.argv[1:])
solutions.read_solutions()
solutions.print()