#!/usr/bin/env python
import glob
import sys
import os
import argparse
import multiprocessing
import subprocess
import codecs
import time
import webbrowser
import shlex

root = os.getcwd()
cpus = multiprocessing.cpu_count()

parser = argparse.ArgumentParser(description="Run ENSO metrics tests",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-H", "--html", action="store_true",
                    help="create and show html result page")
parser.add_argument("-p", "--package", action="store_true",
                    help="package test results")
parser.add_argument(
    "-c",
    "--coverage",
    action="store_true",
    help="run coverage (not implemented)")
parser.add_argument(
    "-v",
    "--verbosity",
    default=1,
    choices=[
        0,
        1,
        2],
    type=int,
    help="verbosity output level")
parser.add_argument(
    "-n",
    "--cpus",
    default=cpus,
    type=int,
    help="number of cpus to use")
parser.add_argument(
    "-f",
    "--failed-only",
    action="store_true",
    default=False,
    help="runs only tests that failed last time and are in the list you provide")
parser.add_argument(
    "-A", "--attributes",
    default=[],
    action="append",
    help="attribute-based runs")
parser.add_argument("tests", nargs="*", help="tests to run")

args = parser.parse_args()


def abspath(path, name, prefix):
    import shutil
    full_path = os.path.abspath(os.path.join(os.getcwd(), "..", path))
    if not os.path.exists(name):
        os.makedirs(name)
    new = os.path.join(name, prefix + "_" + os.path.basename(full_path))
    try:
        shutil.copy(full_path, new)
    except Exception:
        pass
    return new


def findDiffFiles(log):
    i = -1
    file1 = ""
    file2 = ""
    diff = ""
    N = len(log)
    while log[i].find("Source file") == -1 and i > -N:
        i -= 1
    if i > -N:
        file1 = log[i - 1].split()[-1]
        for j in range(i, N):
            if log[j].find("New best!") > -1:
                if log[j].find("Comparing") > -1:
                    file2 = log[j].split()[2]
                else:
                    k = j - 1
                    while log[k].find("Comparing") == -1 and k > -N:
                        k -= 1
                    try:
                        file2 = log[k].split()[2]
                    except Exception:
                        file2 = log[k].split()[1][:-1] + log[j].split()[0]
                        print("+++++++++++++++++++++++++", file2)
            if log[j].find("Saving image diff") > -1:
                diff = log[j].split()[-1]
    return file1, file2, diff


def run_command(command, join_stderr=True):
    if isinstance(command, str):
        command = shlex.split(command)
    if args.verbosity > 0:
        print("Executing %s in %s" % (" ".join(command), os.getcwd()))
    if join_stderr:
        stderr = subprocess.STDOUT
    else:
        stderr = subprocess.PIPE
    P = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=stderr,
        bufsize=0,
        cwd=os.getcwd())
    out = []
    while P.poll() is None:
        read = P.stdout.readline().rstrip().decode("utf-8", errors="replace")
        out.append(read)
        if args.verbosity > 1 and len(read) != 0:
            print(read)
    return P, out


def run_pytest(test_name):
    opts = []
    if args.coverage:
        opts += ["--cov"]
    command = [sys.executable, "-m", "pytest"] + opts + ["-s", test_name]
    start = time.time()
    P, out = run_command(command)
    end = time.time()
    return {test_name: {"result": P.poll(), "log": out, "times": {
        "start": start, "end": end}}}


sys.path.append(
    os.path.join(
        os.path.dirname(
            os.path.abspath(__file__)),
        "tests"))
if len(args.tests) == 0:
    names = glob.glob("tests/test_*.py")
else:
    names = list(args.tests)

if args.failed_only and os.path.exists(os.path.join("tests", ".last_failure")):
    with open(os.path.join("tests", ".last_failure")) as f:
        failed = set(eval(f.read().strip()))
    names = [fnm for fnm in names if fnm in failed]

if args.verbosity > 1:
    print("Names:", names)

if len(names) == 0:
    print("No tests to run")
    sys.exit(0)

p = multiprocessing.Pool(args.cpus)
try:
    outs = p.map_async(run_pytest, names).get(3600)
except KeyboardInterrupt:
    sys.exit(1)
results = {}
failed = []
for d in outs:
    results.update(d)
    nm = next(iter(d))
    if d[nm]["result"] != 0:
        failed.append(nm)
with open(os.path.join("tests", ".last_failure"), "w") as f:
    f.write(repr(failed))

if args.verbosity > 0:
    print("Ran %i tests, %i failed (%.2f%% success)" %
          (len(outs), len(failed), 100. - float(len(failed)) / len(outs) * 100.))
    if len(failed) > 0:
        print("Failed tests:")
        for f in failed:
            print("\t", f)
if args.html or args.package:
    if not os.path.exists("tests_html"):
        os.makedirs("tests_html")
    os.chdir("tests_html")

    fi = open("index.html", "w")
    fi.write("<!DOCTYPE html>\n")
    fi.write("""<html><head><title>ENSO Metrics Test Results %s</title>
    <link rel="stylesheet" type="text/css" href="http://cdn.datatables.net/1.10.13/css/jquery.dataTables.css">
    <script type="text/javascript" src="http://code.jquery.com/jquery-1.12.4.js"></script>
    <script type="text/javascript" charset="utf8"
    src="https://cdn.datatables.net/1.10.13/js/jquery.dataTables.min.js"></script>
    <script>
    $(document).ready( function () {
            $('#table_id').DataTable({
            "order":[[1,'asc'],[0,'asc']],
            "scrollY":"70vh","paging":false,"scrollCollapse":false
            });
                } );
    </script>
    </head>\n""" % time.asctime())
    fi.write("<body><h1>ENSO Metrics Test results: %s</h1>\n" % time.asctime())
    fi.write("<table id='table_id' class='display'>\n")
    fi.write("<thead><tr><th>Test</th><th>Result</th><th>Start Time</th><th>End Time</th><th>Time</th></tr></thead>\n")
    fi.write("<tfoot><tr><th>Test</th><th>Result</th><th>Start Time</th><th>End Time</th><th>Time</th></tr></tfoot>\n")

    for t in sorted(results.keys()):
        result = results[t]
        nm = t.split("/")[-1][:-3]
        fi.write("<tr><td>%s</td>" % nm)
        fe = codecs.open("%s.html" % nm, "w", encoding="utf-8")
        fe.write("<!DOCTYPE html>\n")
        fe.write("<html><head><title>%s</title>" % nm)
        if result["result"] == 0:
            fi.write("<td><a href='%s.html'>OK</a></td>" % nm)
            fe.write("</head><body>")
            fe.write("<a href='index.html'>Back To Results List</a>")
        else:
            fi.write("<td><a href='%s.html'>Fail</a></td>" % nm)
            fe.write("</head><body>")
            fe.write("<a href='index.html'>Back To Results List</a>")
            fe.write("<h1>Failed test: %s on %s</h1>" % (nm, time.asctime()))
        fe.write('<div id="output"><h1>Log</h1><pre>%s</pre></div>' % "\n".join(result["log"]))
        fe.write("<a href='index.html'>Back To Results List</a>")
        fe.write("</body></html>")
        fe.close()
        t_times = result["times"]
        fi.write("<td>%s</td><td>%s</td><td>%s</td></tr>\n" % (
            time.ctime(t_times["start"]), time.ctime(t_times["end"]),
            t_times["end"] - t_times["start"]))

    fi.write("</table></body></html>\n")
    fi.close()
    if args.html:
        webbrowser.open("file://%s/index.html" % os.getcwd())
    os.chdir(root)

if args.package:
    import tarfile
    tnm = "results_%s_%s_%s.tar.bz2" % (
        os.uname()[0], os.uname()[1], time.strftime("%Y-%m-%d_%H:%M"))
    with tarfile.open(tnm, "w:bz2") as t:
        t.add("tests_html")
    if args.verbosity > 0:
        print("Packaged Result Info in:", tnm)

sys.exit(len(failed))
