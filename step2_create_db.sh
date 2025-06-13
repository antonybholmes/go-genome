for f in data/modules/genome/*.sql
do
    out=`echo ${f} | sed -r 's/.sql/.db/'`
    echo ${f} ${out}
    rm ${out}
    cat ${f} | sqlite3 ${out}
    sqlite3 ${out} 'PRAGMA journal_mode=WAL;'
done
