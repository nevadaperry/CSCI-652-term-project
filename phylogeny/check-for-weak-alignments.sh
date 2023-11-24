for f in provided/pw-more/*; do
	echo -n "file: $f pieces: ";
	fgrep -c 'a score=' $f;
done | fgrep -v -e 'pieces: 1'