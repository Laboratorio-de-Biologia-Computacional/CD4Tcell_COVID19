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




f.cfm.table <- function(net, labels, file.name, time=1, replace=list()) {
    cfm <- cellFateMap(net, label.rules=labels, time=time,
                       method="sat.restricted", maxAttractorLength=2) 
    cfm <- apply(cfm,2,as.character)
    cfm <- cfm[  order( cfm[,1], cfm[,2], cfm[,3], cfm[,4] ),]
    write.csv(cfm,file.name, row.names=F)
    cfm <- read.csv(file.name, stringsAsFactors=F)
    for (old in names(replace)) {
        new <- replace[old]
        cfm[cfm[,2]==old,2] <- new
    }
    write.csv(cfm,file.name, row.names=F)
    cfm <- read.csv(file.name, stringsAsFactors=F)
    cfm
}

plot.cfm.alluvial <- function(cfm, file.name=NULL, title='', ignore='Naive') {
    if (!is.null(file.name)) {
        setEPS(); cairo_ps(file.name)
        }
    cfm <- cfm[cfm$initial!=ignore,]
    cfm2d <- cfm %>% group_by(initial, final) %>%
            summarize(freq = n())
    alluvial(cfm2d[,1:2], freq=cfm2d$freq)
    if (!is.null(file.name)) { dev.off() }
    }

plot.node.transitions <- function(cfm, file.name=NULL, title='', normalize=T) {
    cfm.diff <- cfm[cfm$initial!=cfm$final,]
    node.transitions <- merge( summary(as.factor(cfm.diff[cfm.diff$values==0,'genes'])), 
                               summary(as.factor(cfm.diff[cfm.diff$values==1,'genes'])),
                               by=0, all=TRUE)
    rownames(node.transitions) <- node.transitions$Row.names
    node.transitions <- subset(node.transitions, select = -Row.names )
    colnames(node.transitions) <- c('activation','inhibition')
    # format to standarize node names
    node.transitions <- merge( net$genes, node.transitions,
                               by.x=1, by.y=0, all=TRUE)
    node.transitions[is.na(node.transitions)] <- 0
    rownames(node.transitions) <- node.transitions$x
    node.transitions <- subset(node.transitions, select = -x )
    node.transitions <- node.transitions[net$genes,]
    if (normalize) { node.transitions <- node.transitions/summary(as.factor(cfm$genes)) }
    #print(node.transitions)
    # plot plot
    if (!is.null(file.name)) {
        setEPS(); postscript(file.name)
    }
    barplot(t(as.matrix( node.transitions )),
            main=title, xlab='nodes', ylab='% of transitions caused when perturbed',
            col=c('#b3de69','#fb8072'), ylim=c(0,1), las=2, legend=T,
           )
    if (!is.null(file.name)) { dev.off() }
}



f.mut.table <- function(net, lab, file.name) {
    mutants <- perturbNetworkFixedNodes(net, label.rules=lab, returnDataFrame='basinSize')
    mut_socs <- perturbNetworkFixedNodes(net, label.rules=lab, returnDataFrame='basinSize',
                                         genes  = list(c('SOCS1','SOCS2','SOCS3'),c('SOCS1','SOCS2','SOCS3')),
                                         values = list(0,1),
                                         names  = c('SOCSall_0','SOCSall_1'))
    mutants <- merge(mutants, mut_socs, by=0, all=TRUE)
    rownames(mutants) <- mutants$Row.names
    mutants <- mutants[ , -which(names(mutants) %in% c("Row.names"))]
    mutants[mutants == 0] <- NA
    write.csv(mutants,file.name)
    mutants
}

f.mut.plot <- function(mutants, file.name=NULL, title='', normalize=T) {
    if (normalize) {
        mutants <- mutants/mutants
        color <- c('#bebada')
    } else {
        colfunc <- colorRampPalette(c('#bfd3e6', '#810f7c'))
        color <- colfunc(10)
    }
    if (!is.null(file.name)) {
        setEPS(); postscript(file.name)
    }
    mutants
    heatmap(t(as.matrix( mutants )),
            main=title, 
            xlab="Attractor", ylab="Mutant",
            col=color, cexCol=0.75, cexRow=0.75,
            Colv = NA, Rowv = NA, scale="none",
           )
    if (!is.null(file.name)) { dev.off() }

}
