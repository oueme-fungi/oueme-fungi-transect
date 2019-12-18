# Force the first word of the species name to be the genus.
s/Root;(([^;]+;){7})([^;]+);[^; ]+/Root;\1\3;\3/
# Sebacinavermifera is now Serendipta vermifera
s/Sebacina;Sebacina vermifera/Serendipita;Serendipita vermifera/
