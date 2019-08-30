#convert to DNA
/^>/!y/uU/tT/

#fill missing ranks with "unidentified X"
/^>/ { s/Root(;[^;]+)$/&unidentified \\1/
       s/Root(;[^;]+)(;[^;]+)$/&unidentified \\2/;
       s/Root(;[^;]+){{2}}(;[^;]+)$/&unidentified \\2/;
       s/Root(;[^;]+){{3}}(;[^;]+)$/&unidentified \\2/;
       s/Root(;[^;]+){{4}}(;[^;]+)$/&unidentified \\2/;
       s/(unidentified )+/unidentified /g }

#spelling
s/incertae([ _])sedos/incertae\1sedis/g

#Match Unite 
s/Animalia;Metazoa/Metazoa/

#Rosling et al 2011
s/(;Taphrinomycotina incertae sedis){3} -soil clone group 1/;Archaeorhizomycetales;Archaeorhizomycetaceae;Archaeorhizomyces/

#Matheny 2005
s/Cortinariaceae;Inocybe/Inocybaceae;Inocybe/
