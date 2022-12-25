queryGeneDB <- function(table.nm, data.org){
  require('RSQLite')
  
  conv.db <- dbConnect(SQLite(), paste(sqlite.path, data.org, "_genes.sqlite", sep="")); 
  db.map <- dbReadTable(conv.db, table.nm)
  dbDisconnect(conv.db); cleanMem();
  
  return(db.map)
}

.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}
