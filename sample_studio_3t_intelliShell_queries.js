
//db.articles.find({$text: {$search: "dogs"}});


//return related 
db.articles.find({$text: {$search: "dogs"}}, {score: {$meta: "textScore"}}).sort({score:{$meta:"textScore"}}).limit(1000);
