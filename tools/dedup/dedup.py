import sys
infilename = sys.argv[1]
if len(sys.argv) == 3:
  outfilename = sys.argv[2]
else:
  outfilename = infilename + ".cleaned"

print ("in: " + infilename)
print ("out: " + outfilename)

lines_seen = set() # holds lines already seen
outfile = open(outfilename, "w")
for line in open(infilename, "r"):
    if line not in lines_seen: # not a duplicate
        outfile.write(line)
        lines_seen.add(line)
outfile.close()
