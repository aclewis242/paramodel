rm -recurse executable
rm -recurse main.spec
mkdir executable
pyinstaller main.py --distpath executable/dist --workpath executable/build