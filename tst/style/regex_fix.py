#! /usr/bin/env python

import os

subs = [(r"\)\n{", ") {"),
        (r"\){", ") {"),
        (r"(\S),(\w)", r"\1, \2"),
        (r"\t", "  "),
        (r"\n\n\n", "\n\n"),
        #(r"\n\n", "\n"),
        (r"[^#]else\n{", "else {"),
        (r"}\n\s*else", "} else"),
        (r"(\S){", r"\1 {"),
        (r"(\S);(\S)", r"\1; \2"),
        (r"(\w)\s+\[", r"\1["),
        (r"[\s]*$", ""),
        ]


def fix_file(fname, ncols=90):
    import re
    print("Fixing", fname)
    with open(fname, 'r') as f:
        content = f.read()
    out = content
    changed = True
    # apply substitutions
    while changed:
        content = out
        for r, s in subs:
            #print("\t", r, "==>", s)
            out = re.sub(r, s, out)
        changed = content != out
    changed = True
    n = 0
    # fix indentation and line length
    while changed:
        #print("\tFixing indentation and line length")
        content = out
        lines = out.split("\n")
        if lines[-1]:
            lines.append("")
        for i, line in enumerate(lines):
            # handle long lines
            if len(line) > ncols:
                # handle comments
                if re.match(r"\s*//", line):
                    tab = re.match(r"\s+", line)
                    tab = tab.group() if tab else ""
                    comment = "//".join(line.split("//")[1:])
                    comment = re.sub(r"(\S)([*+/-])(\S)", r"\1 \2 \3", comment)
                    words = comment.lstrip().split(" ")
                    buffer = [tab + '//']
                    for word in words:
                        if len(buffer[-1]) + len(word) >= ncols:
                            buffer.append(tab + '// ' + word)
                        else:
                            buffer[-1] += ' ' + word
                    lines[i] = '\n'.join(buffer)
                else:
                    # handle inline comments
                    if "//" in line:
                        tmp = line.split("//")
                        line, comment = tmp[0], "//".join(tmp[1:])
                        tab = re.match("\s+", line)
                        tab = tab.group(0) if tab else ""
                        lines[i] = tab + "//" + comment + "\n" + line
                    # if line asigns a value to a variable
                    elif re.search(r"(?<![<>])=(?!=)", line):
                        loc = re.search(r"(?<![<>])=(?!=)", line).end()
                        prefix = line[:loc]
                        suffix = line[loc:]
                        tab = " " * loc
                        buffer = [prefix]
                        words = [suffix]
                        for delim in " +-*/%":
                            tmp = []
                            for word in words:
                                parts = word.split(delim)
                                for part in parts[:-1]:
                                    tmp.extend([part, delim])
                                tmp.append(parts[-1])
                            words = [i for i in tmp if i]
                        for word in words:
                            if len(buffer[-1]) + len(word) < ncols:
                                buffer[-1] += word
                            else:
                                buffer.append(tab + word.lstrip())
                        buffer[-1] = buffer[-1][:-1]
                        lines[i] = '\n'.join([i.rstrip() for i in buffer])
                    # if line is a function call with multiple arguments
                    elif re.search(r"\w\(.*,.*\)", line):
                        match = re.search(r"([^()]*\w+\(\w*\)\()(.*,.*\))(.*)$", line)
                        if not match:
                            match = re.search(r"([^()]*\w+\()(.*,.*\))(.*)$", line)
                        prefix = match.group(1)
                        suffix = match.group(2) + match.group(3)
                        tab = " " * (len(prefix) - 1)
                        buffer = [prefix]
                        for word in suffix.split(","):
                            if len(buffer[-1]) + len(word) + 1 > ncols:
                                buffer.append(tab + word + ",")
                            else:
                                buffer[-1] += word + ","
                        buffer[-1] = buffer[-1][:-1]
                        lines[i] = '\n'.join([i.rstrip() for i in buffer])
                        if "ROMBERG_TOL" in line:
                            print(line)
                            print("=" * 80)
                            print(lines[i])

                        
        out = "\n".join(lines)
        changed = content != out
        n += 1
        if n > 10:
            for i, line in enumerate(lines):
                if "\n" in line:
                    print("\t\tline", i, "has \\n")
                    print("\t\t", line)
            raise Exception("Too many iterations")
    with open(fname, "w") as f:
        f.write(out)


def fix_dir(dname, ext, **kwargs):
    if isinstance(ext, str):
        ext = [ext]
    for root, dirs, files in os.walk(dname):
        for f in files:
            for e in ext:
                if f.endswith(e):
                    fix_file(os.path.join(root, f), **kwargs)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="?", default=['../../'])
    parser.add_argument("-e", "--extension", nargs="+", default=[".c", ".h"])
    parser.add_argument("-n", "--ncols", type=int, default=90)
    args = parser.parse_args()

    for file in args.files:
        if os.path.isdir(file):
            fix_dir(file, args.extension, ncols=args.ncols)
        else:
            fix_file(file)
