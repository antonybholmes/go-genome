for f in data/modules/genes/*.sql
do
    out=`echo data/modules/genes/${f} | sed -r 's/.sql/.db/'`
    echo ${f} ${out}
    rm ${out}
    cat data/modules/genes/${f} | sqlite3 ${out}
    sqlite3 ${out} 'PRAGMA journal_mode=WAL;'
done
