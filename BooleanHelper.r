get.attr.table <- function(net, labels, file.name=NULL, simplify=F, ...) {
    # calculate attractors
    attr <- getAttractors(net, ...)
    attr.table <- attractorToDataframe(attr, Boolean=TRUE)
    # label by attractor
    attr.labels <- labelAttractors(attr, labels, net_$genes, simplify=simplify)
    attr.table <- merge(x=attr.table, y=attr.labels, by.x = 'attractor', by.y = 0)
    colnames(attr.table)[colnames(attr.table) == "y"] <- "label"
    attr.table$size <- stringr::str_count(attr.table$label,'/')+1
    attr.table <- attr.table[c(c("label", "attractor", "state","size"),net$genes)]
    # sort and save
    if ( !is.null(file.name) ) {
        write.csv(attr.table, file.name, row.names=FALSE)
    }
    attr.table
}

plot.attr.table <- function(attr.table, title='', file.name=NULL, selected.nodes=c()) {
    # select key nodes for plotting
    attr.plot <- subset(attr.table, select = c("label", "attractor", "state", "size", selected.nodes) )
    attr.plot <- attr.plot[with(attr.table, order(label,attractor,state)),]
    
    if (!is.null(file.name)) {
        setEPS(); postscript(file.name)
    }
    attr.plot.steady <- attr.plot[attr.plot$size==1,]
    attr.plot.steady <- subset(attr.plot.steady, select=-c(attractor,state,size))
    attr.plot.steady <- unique(attr.plot.steady)
    heatmap(t(as.matrix( subset(attr.plot.steady, select=-label ))),
            labCol=attr.plot.steady$label,
            main=paste(c(title, "(steady states)"), collapse=' '),
            xlab="Attractor", ylab="Nodes",
            col=c('#fb8072','#b3de69'), cexCol=0.75, cexRow=0.75,
            Colv = NA, Rowv = NA, scale="none",
           )
    if (!is.null(file.name)) { dev.off() }
    
    if (!is.null(file.name)) {
        file.name <- str_replace(file.name, ".eps", "_cycles.eps")
        setEPS(); postscript(file.name)
    }
    attr.plot.cycle <- attr.plot[attr.plot$size>1,]
    if (nrow(attr.plot.cycle)>0) {
    heatmap(t(as.matrix( subset(attr.plot.cycle, select=-c(label,attractor,state,size) ))),
            labCol=attr.plot.cycle$label,
            main=paste(c(title, "(cycles)"), collapse=' '),
            xlab="Attractor", ylab="Nodes",
            col=c('#fb8072','#b3de69'), cexCol=0.75, cexRow=0.75,
            Colv = NA, Rowv = NA, scale="none",
           )
        }
    if (!is.null(file.name)) { dev.off() }
}