#' Our plotting style.
#'
#' @param title The plot title.
#' @param col A vector of colours to use.
#' @param fill A vector of fill colours to use.
#' @param legend_title The legend title.
#' @param legend_position The legend position.
#' @param xlab The title of the x-axis.
#' @param ylab The title of the y-axis.
#' @param font_size The base font size. Default is 11.
#' @return None, add as theme.
AddPlotStyle = function(title=NULL, col=NULL, fill=NULL, legend_title=NULL, legend_position=NULL, xlab=NULL, ylab=NULL, font_size=11) {
  style = list(
    # Basic theme
    ggplot2::theme_light(font_size),
    ggplot2::theme(panel.border = element_blank())
  )
  
  # Title
  if (!is.null(title)) {
    if (nchar(title) > 0) {
      style = c(style, list(ggplot2::ggtitle(title)))
    } else {
      style = c(style, list(ggplot2::theme(title=element_blank())))
    }
  }
  
  # Colour
  if (!is.null(col)) style = c(style, list(ggplot2::scale_colour_manual(values=col)))
  
  # Fill
  if (!is.null(fill)) style = c(style, list(ggplot2::scale_fill_manual(values=fill)))
    
  # Legend title
  if (!is.null(legend_title)) {
    if (nchar(legend_title) > 0) {
      style = c(style, list(labs(color=legend_title, fill=legend_title)))
    } else {
      style = c(style, list(ggplot2::theme(legend.title=element_blank())))
    }
  }
  
  # Legend position
  if (!is.null(legend_position)) style = c(style, list(ggplot2::theme(legend.position=legend_position)))
  
  # Axis labels
  if (!is.null(xlab)) {
    if (nchar(xlab) > 0) {
      style = c(style, list(ggplot2::xlab(xlab)))
    } else {
      style = c(style, list(ggplot2::theme(axis.title.x=element_blank())))
    }
  }
  
  if (!is.null(ylab)) {
    if (nchar(ylab) > 0) {
      style = c(style, list(ggplot2::ylab(ylab)))
    } else {
      style = c(style, list(ggplot2::theme(axis.title.y=element_blank())))
    }
  }
  
  return(style)
}

#' Helper function to generate plot captions.
#'
#' @param plot_names A list of plot names.
#' @param assay_names A list of assay names that should be recognized in the plot names as such.
#' @param split If not NULL, split plot names to generate a X vs Y caption.
#' @param capitalize If TRUE, capitalize the first letter of the caption.
#' @return A list of captions.
GeneratePlotCaptions = function(plot_names, assay_names=NULL, split=NULL, capitalize=TRUE) {
  if (length(plot_names) == 0) return(NULL)

  # Split names into two parts if requested
  if (!is.null(split)) {
    plot_names = strsplit(plot_names, split)
  } else {
    plot_names = as.list(plot_names)
  }
  
  # Match and generate captions
  captions = purrr::map(plot_names, function(nms) {
      
    # Add captions here
    capts = dplyr::case_when(startsWith(nms, "nCount_") ~ "number of counts",
                             startsWith(nms, "nFeature_") ~ "number of features",
                             startsWith(nms, "pCountsTop50_") ~ "percent counts in top50 features",
                             startsWith(nms, "pCountsTop5_") ~ "percent counts in top5 features",
                             startsWith(nms, "pMito_") ~ "percent counts in mitochondrial genes",
                             startsWith(nms, "pRibosomal_") ~ "percent counts in ribosomal genes",
                             startsWith(nms, "pGlobin_") ~ "percent counts in globin genes",
                             startsWith(nms, "pERCC_") ~ "percent counts in ERCC controls",
                             startsWith(nms, "pXIST_") ~ "percent counts in female-determining XIST gene",
                             startsWith(nms, "pChrY_") ~ "percent counts in male-determining chrY genes",
                             startsWith(nms, "S.Score") ~ "S phase score",
                             startsWith(nms, "G2M.Score") ~ "G2M phase score",
                             startsWith(nms, "Phase") ~ "cell cycle phase",
                             startsWith(nms, "varFeatures") ~ "variable features",
                             startsWith(nms, "Expression") ~ "Expression of",
                             .default = NULL
    )
    
    # Add generic caption if necessary
    idx = which(is.na(capts))
    capts[idx] = paste0("metric '", nms[idx], "'")
    
    # Try to find out whether the assay is part of the name (or its parts) and add it to the caption
    if (!is.null(assay_names)) {
        
        # Iterate over the plot names
        assays_to_add = purrr::map_chr(nms, function(n) {
            # Check which assay name is part of the name (or its parts) 
            is_part = purrr::map_lgl(assay_names, function(a) return(endsWith(n, paste0("_", a))))
            
            if (any(is_part)) {
                # If yes, return the assay name
                return(assay_names[is_part])
            } else {
                # If not, return NA
                return(NA)
            }
        })
        
        # If all parts of the name contain the same assay name, add it to the last part else to all parts
        if (length(unique(na.omit(assays_to_add))) == 1) {
            i = length(capts)
            
        } else {
            i = which(!is.na(assays_to_add))
        }
        capts[i] = paste0(capts[i], " for assay '", assays_to_add[i], "'")
    }
    
    return(capts)
  })
  
  # If there are two parts, collapse with " Vs "
  captions = purrr::map(captions, function(x) {
    for(i in seq_along(x)) {
      # Skip first one
      if (i==1) next
      
      # Make first character lower-case unless the first two characters are upper-case (meaning an abbreviation)
      # if (!grepl("^[[:upper:]]{2,}", x[i])) {
      #   s = unlist(strsplit(x[i], ""))
      #   s[1] = tolower(s[1])
      #   x[i] = paste(s, collapse="")
      # }
    }
    
    # Collapse
    x = paste(x, collapse=" vs ")
    
    # If requested, make first letter upper-case
    if (capitalize) {
      s = unlist(strsplit(x, ""))
      s[1] = toupper(s[1])
      x = paste(s, collapse="")
    }
    
    return(x)
  }) %>% unlist()
  
  return(captions)
}

#' Plots barcode metadata for QC. Numeric columns will plotted as violin plots and non-numeric
#' columns will be plotted as bar plots.
#'
#' @param sc Seurat v5 object.
#' @param qc Barcode metadata columns to plot.
#' @param filter A nested list where the first level is the barcode metadata column and the second levels 
#' contains filters per dataset. Filters for numeric columns must numeric vectors with min and max. Filter
#' for character/factor columns must be character vectors with the values that should be kept. 
#' level contains the filter values
#' @param assay The assay of the barcodes. If NULL, defaults to default assay of the Seurat object.
#' @param log10 If set to TRUE, show the y-axis in log10. Can also be a regular expression to apply only to specific columns. Only numeric columns will be affected.
#' @return A list of ggplot2 objects.
PlotBarcodeQC = function(sc, qc, filter=NULL, assay=NULL, log10=FALSE) {
  # Get assay and barcodes
  if (is.null(assay)) assay = SeuratObject::DefaultAssay(sc)
  bcs = SeuratObject::Cells(sc[[assay]])
  barcode_metadata = sc[[]][bcs, ]
  
  # Determine QC type (numeric or non-numeric)
  is_numeric = purrr::map_lgl(qc, function(q) return(is.numeric(barcode_metadata[, q, drop=TRUE])))
  
  # Get filter thresholds per QC metrics (numeric)
  qc_thresholds_numeric = purrr::map(qc[is_numeric], function(f) {
    tresh = purrr::map_dfr(names(filter), function(n) {
      tr = data.frame(qc_feature=character(), ident=character(), 
                      threshold=character(), value=numeric())
      
      if (f %in% names(filter[[n]])) {
        tr = data.frame(qc_feature=f, 
                        ident=n, 
                        min=filter[[n]][[f]][1], 
                        max=filter[[n]][[f]][2])  %>% 
          tidyr::pivot_longer(c(min, max), names_to="threshold", values_to="value")
      }
      
      tr$ident = factor(tr$ident, levels=orig_idents)
      return(tr)
    })
    return(tresh)
  })
  names(qc_thresholds_numeric) = qc[is_numeric]
  
  qc_thresholds_other = purrr::map(qc[!is_numeric], function(f) {
    tresh = purrr::map(names(filter), function(n) {
      return(filter[[n]][[f]])
    })
    names(tresh) = names(filter)
    return(tresh)
  })
  names(qc_thresholds_other) = qc[!is_numeric]
  
  # Plot numeric QC
  plist_numeric = Seurat::VlnPlot(sc, features=qc[is_numeric], combine=FALSE, pt.size=0, layer="counts", assay=assay)
  names(plist_numeric) = qc[is_numeric]
  
  plist_numeric = purrr::map(names(plist_numeric), function(n) {
    # Add style
    p = plist_numeric[[n]] + 
      AddPlotStyle(legend_position="none", xlab="", fill=ScColours(sc, "orig.ident")) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    # Add filter thresholds
    qc_threshold_segments = purrr::pmap(qc_thresholds_numeric[[n]], function(qc_feature, ident, threshold, value) {
      return(annotate(geom="segment", x=as.integer(ident)-0.5, xend=as.integer(ident)+0.5, y=value, yend=value, linetype=ifelse(threshold=="max", 6, 2), colour="firebrick"))
    })
    p = p + qc_threshold_segments
    
    
    # Add log10 transformation if requested
    if (!is.null(log10)) {
      if ( (is.logical(log10) & log10==TRUE) | (is.character(log10) & grepl(pattern=log10, x=n)) ) {
          p = p + scale_y_continuous(trans="log10")
      }
    }
    return(p)
  })
  names(plist_numeric) = qc[is_numeric]
  
  # Plot non-numeric QC
  plist_other = purrr::map(qc[!is_numeric], function(n) {
    # Collect plot_data
    plot_data = barcode_metadata %>% 
      dplyr::select(x=orig.ident, y=!!sym(n)) %>%
      dplyr::count(x, y) %>%
      dplyr::group_by(x) %>%
      dplyr::mutate(perc=n/sum(n)*100)
    plot_data$status = "filter"
    for(i in 1:nrow(plot_data)){
      allowed_values = qc_thresholds_other[[n]]
      allowed_values = allowed_values[[plot_data$x[i]]]
      plot_data$status[i] = ifelse(plot_data$y[i] %in% allowed_values, "keep", "filter")
    }
    plot_data$status = factor(plot_data$status, levels=c("keep", "filter"))
    
    # Make plot
    p = ggplot(plot_data, aes(x=x, y=perc, fill=y, alpha=status, col=status)) +
      geom_bar(stat="identity", position="stack") +
      scale_x_discrete("Identity") +
      scale_color_manual(values=c("keep"="black", "filter"="grey")) +
      scale_alpha_manual(values=c("keep"=1, "filter"=0.2)) +
      AddPlotStyle(fill=ScColours(sc, "orig.ident")) +
      theme(axis.title.y=element_blank())
    
    return(p)
  })
  names(plist_other) = qc[!is_numeric]
  
  # Adjust order
  plist = c(plist_numeric, plist_other)
  plist = plist[qc]
  
  return(plist)
}

#' Plots two barcode metadata columns for QC. Supports only numeric columns.
#'
#' @param sc Seurat v5 object.
#' @param qc Pairs of barcode metadata columns to plot.
#' @param filter A nested list where the first level is the barcode metadata column and the second levels 
#' contains filters per dataset. Filters for numeric columns must numeric vectors with min and max. Filter
#' @param assay The assay of the barcodes. If NULL, defaults to default assay of the Seurat object.
#' for character/factor columns must be character vectors with the values that should be kept.
#' @return A list of ggplot2 objects.
PlotBarcodeQCCor = function(sc, qc, filter=NULL, assay=NULL) {
  # Get assay and barcodes
  if (is.null(assay)) assay = SeuratObject::DefaultAssay(sc)
  bcs = SeuratObject::Cells(sc[[assay]])
  barcode_metadata = sc[[]][bcs, ]
  
  # Check QC type (only numeric allowed)
  qc_cols = purrr::flatten(qc) %>%
    unlist() %>%
    unique()
  f = purrr::map_lgl(qc_cols, function(c) return(is.numeric(barcode_metadata[, c, drop=TRUE])))
  assertthat::assert_that(all(f),
                          msg=FormatString("Barcode metadata columns {qc_cols[!f]*} are not numeric! Function PlotBarcodeQCCor can only plot numeric data."))
  
  # Get filter thresholds per QC metrics (numeric)
  qc_thresholds = purrr::map(qc_cols, function(f) {
    tresh = purrr::map_dfr(names(filter), function(n) {
      tr = data.frame(qc_feature=character(), ident=character(), 
                      threshold=character(), value=numeric())
      
      if (f %in% names(filter[[n]])) {
        tr = data.frame(qc_feature=f, 
                        ident=n, 
                        min=filter[[n]][[f]][1], 
                        max=filter[[n]][[f]][2])  %>% 
          tidyr::pivot_longer(c(min, max), names_to="threshold", values_to="value")
      }
      
      tr$ident = factor(tr$ident, levels=orig_idents)
      return(tr)
    })
    return(tresh)
  })
  names(qc_thresholds) = qc_cols
  
  plist = purrr::map(qc, function(c) {
    f1 = c[1]
    f2 = c[2]
    
    # Plot QC feature f1 vs f2
    p = Seurat::FeatureScatter(sc, cells=SeuratObject::Cells(sc[[assay]]), feature1=f1, feature2=f2, shuffle=TRUE, seed=getOption("random_seed"))
    p = p + AddPlotStyle(col=ScColours(sc, "orig.ident"))
    
    # Add filter thresholds for f1
    qc_threshold_segments = purrr::pmap(qc_thresholds[[f1]], function(qc_feature, ident, threshold, value) {
      return(geom_vline(xintercept=value, linetype=ifelse(threshold=="max", 6, 2), colour="firebrick"))
    })
    p = p + qc_threshold_segments
    
    # Add filter thresholds for f2
    qc_threshold_segments = purrr::pmap(qc_thresholds[[f2]], function(qc_feature, ident, threshold, value) {
      return(geom_hline(yintercept=value, linetype=ifelse(threshold=="max", 6, 2), colour="firebrick"))
    })
    p = p + qc_threshold_segments
    
    return(p)
  })
  names(plist) = purrr::map(qc, paste, collapse="__") %>% unlist()
  
  return(plist)
}

#' Plots the variable features for each layer (dataset).
#'
#' @param sc Seurat v5 object.
#' @param method Method used to find variable features. Can be: 'vst' (Seurat standard), 'sct' (SCTransform) or 'scran' (Scran).
#' @param assay Assay. If NULL, defaults to default assay of the Seurat object.
#' @param top The top genes that should be labeled.
#' @return A list of ggplot2 objects.
PlotVariableFeatures = function(sc, method, assay=NULL, top=10) {

  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Checks
  valid_methods = c("vst", "sct", "scran")
  assertthat::assert_that(method %in% valid_methods,
                          msg=FormatString("Method is {method} but must be one of: {valid_methods*}."))
  
  layers = SeuratObject::Layers(sc[[assay]], "data")
  assertthat::assert_that(length(layers) > 0,
                          msg=FormatString("Could not find normalized data for assay {assay}."))
  
  # Make plots per layer (dataset)
  orig_idents =  levels(sc$orig.ident)
  
  plist = purrr::map(orig_idents, function(n) {
      
    # Collect information about highly variable genes
    if (method == "sct") {
      assertthat::assert_that(.hasSlot(sc[[assay]], "SCTModel.list"),
                              msg=FormatString("No variable feature information in slot SCTModel.list available for {assay}."))
      hvf_info = Seurat::SCTResults(sc[[assay]], slot="feature.attributes", model=n)
      hvf_info = hvf_info[, c("gmean", "variance", "residual_variance")]
      hvf_info$variable = rownames(hvf_info) %in% Seurat::VariableFeatures(sc[[assay]])
      hvf_info$rank = match(rownames(hvf_info), Seurat::VariableFeatures(sc[[assay]]))
    } else if (method == "vst") {
      hvf_info = SeuratObject::HVFInfo(sc[[assay]],
                                       method="vst",
                                       layer=paste("counts", n, sep="."),
                                       status=TRUE)
      if (is.null(hvf_info)) {
        hvf_info = SeuratObject::HVFInfo(sc[[assay]],
                                 method="vst",
                                 layer=paste("counts", n, sep="."),
                                 status=TRUE)
      }
      assertthat::assert_that(!is.null(hvf_info),
                              msg=FormatString("No variable feature information available for {assay}."))
    } else if (method == "scran") {
      hvf_info = SeuratObject::HVFInfo(sc[[assay]],
                                       method="scran",
                                       layer=paste("counts", n, sep="."),
                                       status=TRUE)
      if (is.null(hvf_info)) {
        hvf_info = SeuratObject::HVFInfo(sc[[assay]],
                                 method="scran",
                                 layer=paste("counts", n, sep="."),
                                 status=TRUE)
      }
      assertthat::assert_that(!is.null(hvf_info),
                              msg=FormatString("No variable feature information available for {assay}."))
      
    }
    
    # Get top genes
    top_genes = which(hvf_info$rank<=top)
    top_genes = rownames(hvf_info)[top_genes]
    
    # Define columns to plot
    if (method == "scran") {
      hvf_info = hvf_info[, c("mean", "total", "variable", "rank")]
      xlab = "Average Expression"
      ylab = "Dispersion"
    } else if (method == "vst") {
      
      hvf_info = hvf_info[, c("mean", "variance.standardized", "variable", "rank")]
      xlab = "Average Expression"
      ylab = "Standardized Variance"
    } else if (method == "sct") {
      hvf_info = hvf_info[, c("gmean", "residual_variance", "variable", "rank")]
      xlab = "Geometric Mean of Expression"
      ylab = "Residual Variance"
    }
    
    xvar = colnames(hvf_info)[1]
    yvar = colnames(hvf_info)[2]
    
    # Make plot
    p = ggplot(hvf_info, aes(x=!!sym(xvar), y=!!sym(yvar), col=variable)) +
      geom_point() +
      scale_x_log10(xlab) +
      scale_y_continuous(ylab) +
      scale_color_manual("Variable", values=c("FALSE"="black", "TRUE"="red")) +
      AddPlotStyle() + 
      theme(legend.position="none", legend.background=element_rect(fill=alpha("white", 0.0)))
    p = Seurat::LabelPoints(plot=p, points=top_genes, repel=TRUE, xnudge=0, ynudge=0)
    return(p)
  })
  names(plist) = paste("varFeatures", orig_idents, sep="_")
  
  return(plist)
}

#' Plots the relative log expression features for each layer (dataset).
#'
#' @param sc Seurat v5 object.
#' @param assay Assay. If NULL, defaults to default assay of the Seurat object.
#' @param layer Type of layer. Can be counts or data but also a specific layer.
#' @param nbarcodes Number of barcodes to plot.
#' @param is_log Is the data already in log? If not, will be logged.
#' @return A ggplot2 object.
PlotRLE = function(sc, assay=NULL, layer="counts", nbarcodes=500, is_log=FALSE) { 
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Checks
  layers = SeuratObject::Layers(sc[[assay]], layer)
  assertthat::assert_that(length(layers) > 0,
                          msg=FormatString("Could not find data for layer {layer} of assay {assay}."))
  
  if (!is(sc[[assay]], "SCTAssay")) {
    # Standard assays
    
    # Sample
    data = purrr::map(layers, function(l) {
      dt = SeuratObject::LayerData(sc, assay=assay, layer=l)
      set.seed(getOption("random_seed"))
      idx = 1:ncol(dt)
      smp = sample(idx, min(length(idx), nbarcodes))
      return(dt[, smp])
    })
    
    # Get identity
    orig_idents = purrr::map(seq_along(layers), function(i) {
      return(rep(layers[i], ncol(data[[i]])))
    }) %>% unlist()
    orig_idents = gsub(pattern=paste0(layer,"."), replacement="", x=orig_idents)
    orig_idents = factor(orig_idents, levels=levels(sc$orig.ident))
    
    # Merge
    data = purrr::reduce(data, cbind)
  } else {
    # SCT assay
    
    # Get data
    data = SeuratObject::LayerData(sc, assay=assay, layer=layer)
    
    # Sample
    barcode_metadata = sc[[]]
    barcode_metadata = barcode_metadata[, ]
    orig_idents = split(colnames(data), barcode_metadata$orig.ident)
    orig_idents = purrr::map(orig_idents, function(bcs) {
      set.seed(getOption("random_seed"))
      bcs = sample(bcs, min(length(bcs), nbarcodes))
      return(bcs)
    })
    
    # Subset data
    data = data[, orig_idents %>% purrr::flatten_chr()]
    data = as(data, "dgCMatrix")
    
    # Get identity
    orig_idents = purrr::map(names(orig_idents), function(n) {
      return(rep(n, length(orig_idents[[n]])))
    }) %>% purrr::flatten_chr()
    orig_idents = factor(orig_idents, levels=levels(sc$orig.ident))
  }
    
  if (!is_log) {
    data = log1p(data + 1)
  }
  
  # Calculate feature medians across all barcodes
  row_medians = CalculateMedians(matrix=data, margin=1, chunk_size=NULL)
  
  # Subtract gene median from gene count
  data = data - row_medians
  
  # Calculate RLE stats per cell
  rle_stats = CalculateBoxplotStats(data, margin=2, chunk_size=NULL)
  rle_stats$orig.ident = orig_idents
  rle_stats = rle_stats %>% 
    dplyr::arrange(orig.ident)
  rle_stats$x = 1:nrow(rle_stats)
  
  
  # Actual plotting
  p = ggplot(rle_stats) +
    geom_segment(aes(x=x, xend=x, y=lower_whisker , yend=upper_whisker, col=orig.ident)) +
    geom_segment(aes(x=x, xend=x, y=q25-0.01 , yend=q75+0.01), colour="grey20") +
    geom_point(aes(x=x, y=q50), shape=1) +
    AddPlotStyle(xlab="Cells", ylab="Relative log expression", col=ScColours(sc, "orig.ident"), legend_title="Dataset") + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  return(p)
}

#' Wrapper for spatial dim plots. Takes as input a Seurat v5 object, one or more image names and other parameters to 
#' be passed on to SpatialDimPlot (sequencing-based) or ImageDimPlot (image-based).
#'
#' @param sc Seurat v5 object.
#' @param images One or more images. If NULL, will use all images in Seurat object.
#' @param assay Get all images with this default assay.
#' @return A list of ggplot2 objects.
DimPlotSpatial = function(sc, assay=NULL, images=NULL, ...) {
  # If NULL, use the default assay of the Seurat object
  if (is.null(assay)) assay = SeuratObject::DefaultAssay(sc)
  
  # If images is NULL, get all images with default assay 'image_assay'
  if (is.null(images)) {
    images = SeuratObject::Images(sc, assay=assay)
  }

  # Plot
  plist = purrr::map(images, function(i) {
    image_type = class(sc@images[[i]])
    if (grepl("Visium", image_type)) {
      # Visium image plot
      return(Seurat::SpatialDimPlot(sc, images=i, ...))
    } else {
      # Xenium FOV plot
      return(Seurat::ImageDimPlot(sc, fov=i, ...))
    }
  })
  names(plist) = images
  return(plist)
}

#' Wrapper for spatial feature plots. Takes as input a Seurat v5 object, one or more image names and other parameters to 
#' be passed on to SpatialFeaturePlot (sequencing-based) or ImageFeaturePlot (image-based).
#'
#' @param sc Seurat v5 object.
#' @param images One or more images. If NULL, will use all images in Seurat object.
#' @param assay Get all images with this default assay.
#' @return A list of ggplot2 objects.
FeaturePlotSpatial = function(sc, assay=NULL, images=NULL, ...) {
  # If NULL, use the default assay of the Seurat object
  if (is.null(assay)) assay = SeuratObject::DefaultAssay(sc)
  
  # If images is NULL, get all images with default assay 'image_assay'
  if (is.null(images)) {
    images = SeuratObject::Images(sc, assay=assay)
  }
  
  # Plot
  plist = purrr::map(images, function(i) {
    image_type = class(sc@images[[i]])
    if (grepl("Visium", image_type)) {
      # Visium image plot
      return(Seurat::SpatialFeaturePlot(sc, images=i, ...))
    } else {
      # Xenium FOV plot
      return(Seurat::ImageFeaturePlot(sc, fov=i, ...))
    }
  })
  names(plist) = images
  return(plist)
}

#' Transform a matrix cells (rows) x htos (cols) into a format that can be understood by feature_grid: cell class, name hto1, value hto1, name hto2, value hto2
#' 
#' @param x: A matrix cells (rows) x htos (cols).
#' @param cell_classification A vector of cell classifications.
DfAllColumnCombinations = function(x, cell_classification) {
  out = combn(x, 2, simplify=FALSE)
  out = lapply(out, function(o) {
    return(data.frame(cell_classification=unname(cell_classification[rownames(o)]), name1=colnames(o)[1], value1=o[, 1], name2=colnames(o)[2], value2=o[, 2]))
  }) %>% dplyr::bind_rows()
  
  # Define plot order so that the two levels of interest are always on top, then negatives, doublets, 
  #   and finally all other samples
  out$order = 0
  out[out$cell_classification==out$name1, "order"] = 3
  out[out$cell_classification==out$name2, "order"] = 3
  out[out$cell_classification=="Negative", "order"] = 2
  out[out$cell_classification=="Doublet", "order"] = 1
  out = out %>% dplyr::group_by(name1, name2) %>% dplyr::arrange(order)
  out$order = NULL # remove column again -> only needed to order data points
  
  return(out)
}

################################################################################
# IMPORTANT: The following functions are part of the pull request that implements
# reading of spaceranger segmentations for Seurat: https://github.com/satijalab/seurat/pull/10028
################################################################################

# Once they are part of an official Seurat release, this code needs to be removed.
assertthat::assert_that(exists("SpatialFeaturePlot1") | !"plot_segmentations" %in% names(formals(Seurat::SpatialPlot)),
                        msg="Seurat::SpatialFeaturePlot now supports plotting segmentations. Please remove all function definitions below in R/functions_plotting.R!")
assertthat::assert_that(exists("SpatialDimPlot1") | !"plot_segmentations" %in% names(formals(Seurat::SpatialPlot)),
                        msg="Seurat::SpatialDimPlot now supports plotting segmentations. Please remove all function definitions below in R/functions_plotting.R!")

SpatialFeaturePlot1 <- function(
        object,
        features,
        images = NULL,
        crop = TRUE,
        slot = 'data',
        keep.scale = "feature",
        min.cutoff = NA,
        max.cutoff = NA,
        ncol = NULL,
        combine = TRUE,
        pt.size.factor = 1.6,
        alpha = c(1, 1),
        image.alpha = 1,
        image.scale = "lowres",
        shape = 21,
        stroke = NA,
        interactive = FALSE,
        information = NULL,
        plot_segmentations = FALSE
) {
    return(SpatialPlot(
        object = object,
        features = features,
        images = images,
        crop = crop,
        slot = slot,
        keep.scale = keep.scale,
        min.cutoff = min.cutoff,
        max.cutoff = max.cutoff,
        ncol = ncol,
        combine = combine,
        pt.size.factor = pt.size.factor,
        alpha = alpha,
        image.alpha = image.alpha,
        image.scale = image.scale,
        shape = shape,
        stroke = stroke,
        interactive = interactive,
        information = information,
        plot_segmentations = plot_segmentations
    ))
}
assignInNamespace("SpatialFeaturePlot",SpatialFeaturePlot1,ns="Seurat")
SpatialFeaturePlot = SpatialFeaturePlot1

SpatialDimPlot1 <- function(
        object,
        group.by = NULL,
        images = NULL,
        cols = NULL,
        crop = TRUE,
        cells.highlight = NULL,
        cols.highlight = c('#DE2D26', 'grey50'),
        facet.highlight = FALSE,
        label = FALSE,
        label.size = 7,
        label.color = 'white',
        repel = FALSE,
        ncol = NULL,
        combine = TRUE,
        pt.size.factor = 1.6,
        alpha = c(1, 1),
        image.alpha = 1,
        image.scale = "lowres",
        shape = 21,
        stroke = NA,
        label.box = TRUE,
        interactive = FALSE,
        information = NULL,
        plot_segmentations = FALSE
) {
    return(SpatialPlot(
        object = object,
        group.by = group.by,
        images = images,
        cols = cols,
        crop = crop,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        facet.highlight = facet.highlight,
        label = label,
        label.size = label.size,
        label.color = label.color,
        repel = repel,
        ncol = ncol,
        combine = combine,
        pt.size.factor = pt.size.factor,
        alpha = alpha,
        image.alpha = image.alpha,
        image.scale = image.scale,
        shape = shape,
        stroke = stroke,
        label.box = label.box,
        interactive = interactive,
        information = information,
        plot_segmentations = plot_segmentations
    ))
}
assignInNamespace("SpatialDimPlot",SpatialDimPlot1,ns="Seurat")
SpatialDimPlot = SpatialDimPlot1

SpatialPlot1 <- function(
        object,
        group.by = NULL,
        features = NULL,
        images = NULL,
        cols = NULL,
        image.alpha = 1,
        image.scale = "lowres",
        crop = TRUE,
        slot = 'data',
        keep.scale = "feature",
        min.cutoff = NA,
        max.cutoff = NA,
        cells.highlight = NULL,
        cols.highlight = c('#DE2D26', 'grey50'),
        facet.highlight = FALSE,
        label = FALSE,
        label.size = 5,
        label.color = 'white',
        label.box = TRUE,
        repel = FALSE,
        ncol = NULL,
        combine = TRUE,
        pt.size.factor = 1.6,
        alpha = c(1, 1),
        shape = 21,
        stroke = NA,
        interactive = FALSE,
        do.identify = FALSE,
        identify.ident = NULL,
        do.hover = FALSE,
        information = NULL,
        plot_segmentations = FALSE
) {
    if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
        warning(
            "'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity",
            call. = FALSE,
            immediate. = TRUE
        )
        interactive <- TRUE
    }
    if (!is.null(x = group.by) & !is.null(x = features)) {
        stop("Please specific either group.by or features, not both.")
    }
    images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
    if (length(x = images) == 0) {
        images <- Images(object = object)
    }
    if (length(x = images) < 1) {
        stop("Could not find any spatial image information")
    }
    
    # Check keep.scale param for valid entries
    if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", "all"))) {
        stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
    }
    
    cells <- unique(CellsByImage(object, images = images, unlist = TRUE))
    if (is.null(x = features)) {
        if (interactive) {
            return(ISpatialDimPlot(
                object = object,
                image = images[1],
                image.scale = image.scale,
                group.by = group.by,
                alpha = alpha
            ))
        }
        group.by <- group.by %||% 'ident'
        object[['ident']] <- Idents(object = object)
        data <- object[[group.by]]
        data <- data[cells,,drop=F]
        for (group in group.by) {
            if (!is.factor(x = data[, group])) {
                data[, group] <- factor(x = data[, group])
            }
        }
    } else {
        if (interactive) {
            return(ISpatialFeaturePlot(
                object = object,
                feature = features[1],
                image = images[1],
                image.scale = image.scale,
                slot = slot,
                alpha = alpha
            ))
        }
        data <- FetchData(
            object = object,
            vars = features,
            cells = cells,
            layer = slot,
            clean = FALSE
        )
        features <- colnames(x = data)
        # Determine cutoffs
        min.cutoff <- mapply(
            FUN = function(cutoff, feature) {
                return(ifelse(
                    test = is.na(x = cutoff),
                    yes = min(data[, feature]),
                    no = cutoff
                ))
            },
            cutoff = min.cutoff,
            feature = features
        )
        max.cutoff <- mapply(
            FUN = function(cutoff, feature) {
                return(ifelse(
                    test = is.na(x = cutoff),
                    yes = max(data[, feature]),
                    no = cutoff
                ))
            },
            cutoff = max.cutoff,
            feature = features
        )
        check.lengths <- unique(x = vapply(
            X = list(features, min.cutoff, max.cutoff),
            FUN = length,
            FUN.VALUE = numeric(length = 1)
        ))
        if (length(x = check.lengths) != 1) {
            stop("There must be the same number of minimum and maximum cuttoffs as there are features")
        }
        # Apply cutoffs
        data <- sapply(
            X = 1:ncol(x = data),
            FUN = function(index) {
                data.feature <- as.vector(x = data[, index])
                min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
                max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
                data.feature[data.feature < min.use] <- min.use
                data.feature[data.feature > max.use] <- max.use
                return(data.feature)
            }
        )
        colnames(x = data) <- features
        rownames(x = data) <- cells
    }
    features <- colnames(x = data)
    colnames(x = data) <- features
    rownames(x = data) <- cells
    facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && is.list(x = cells.highlight))
    if (do.hover) {
        if (length(x = images) > 1) {
            images <- images[1]
            warning(
                "'do.hover' requires only one image, using image ",
                images,
                call. = FALSE,
                immediate. = TRUE
            )
        }
        if (length(x = features) > 1) {
            features <- features[1]
            type <- ifelse(test = is.null(x = group.by), yes = 'feature', no = 'grouping')
            warning(
                "'do.hover' requires only one ",
                type,
                ", using ",
                features,
                call. = FALSE,
                immediate. = TRUE
            )
        }
        if (facet.highlight) {
            warning(
                "'do.hover' requires no faceting highlighted cells",
                call. = FALSE,
                immediate. = TRUE
            )
            facet.highlight <- FALSE
        }
    }
    if (facet.highlight) {
        if (length(x = images) > 1) {
            images <- images[1]
            warning(
                "Faceting the highlight only works with a single image, using image ",
                images,
                call. = FALSE,
                immediate. = TRUE
            )
        }
        ncols <- length(x = cells.highlight)
    } else {
        ncols <- length(x = images)
    }
    plots <- vector(
        mode = "list",
        length = length(x = features) * ncols
    )
    
    # Get max across all features
    if (!(is.null(x = keep.scale)) && keep.scale == "all") {
        max.feature.value <- max(apply(data, 2, function(x) max(x, na.rm = TRUE)))
    }
    
    for (i in 1:ncols) {
        plot.idx <- i
        image.idx <- ifelse(test = facet.highlight, yes = 1, no = i)
        image.use <- object[[images[[image.idx]]]]
        #Extract image information
        
        coordinates <- GetTissueCoordinates(
            object = image.use,
            scale = image.scale
        )
        #CRITICAL STEP: if the rownames do not match the cell ids, then dataframe is not created properly
        rownames(coordinates) <- coordinates$cell
        highlight.use <- if (facet.highlight) {
            cells.highlight[i]
        } else {
            cells.highlight
        }
        for (j in 1:length(x = features)) {
            cols.unset <- is.factor(x = data[, features[j]]) && is.null(x = cols)
            if (cols.unset) {
                cols <- scales::hue_pal()(n = length(x = levels(x = data[, features[j]])))
                names(x = cols) <- levels(x = data[, features[j]])
            }
            
            # Get feature max for individual feature
            if (!(is.null(x = keep.scale)) && keep.scale == "feature" && !inherits(x = data[, features[j]], what = "factor") ) {
                max.feature.value <- max(data[, features[j]])
            }
            
            #WARNING: The dataframe creation step takes a long time
            #Has been shown to take upwards of 5 minutes
            #Positional indexing
            common_cells <- intersect(rownames(coordinates), rownames(data))
            coord_idx <- match(common_cells, rownames(coordinates))
            data_idx <- match(common_cells, rownames(data))
            
            dataframe <- cbind(
                coordinates[coord_idx, ],
                data[data_idx, features[j], drop = FALSE]
            )
            
            #Check if object contains a sf slot (attached via Load10X_Spatial)
            #If so, use sf-geometry based rendering 
            use_geom_sf <- (inherits(image.use, "VisiumV2") &&
                                !is.null(image.use@boundaries$segmentation) &&
                                "sf.data" %in% slotNames(image.use@boundaries$segmentation))
            
            #WARNING: The dataframe creation step takes a long time
            plot <- SingleSpatialPlot(
                data = cbind(
                    coordinates,
                    data[rownames(x = coordinates), features[j], drop = FALSE]
                ),
                image = image.use,
                image.scale = image.scale,
                image.alpha = image.alpha,
                col.by = features[j],
                cols = cols,
                alpha.by = if (is.null(x = group.by)) {
                    features[j]
                } else {
                    NULL
                },
                pt.alpha = if (!is.null(x = group.by)) {
                    alpha[j]
                } else {
                    NULL
                },
                geom = if (inherits(x = image.use, what = "STARmap")) {
                    'poly'
                } else if (use_geom_sf) {
                    # Use sf for both segmentations and centroids when sf data is available
                    "sf"
                } else {
                    "spatial"
                },
                cells.highlight = highlight.use,
                cols.highlight = cols.highlight,
                pt.size.factor = pt.size.factor,
                shape = shape,
                stroke = stroke,
                crop = crop,
                plot_segmentations = plot_segmentations
            )
            if (is.null(x = group.by)) {
                plot <- plot +
                    scale_fill_gradientn(
                        name = features[j],
                        colours = SpatialColors(n = 100)
                    ) +
                    theme(legend.position = 'top') +
                    scale_alpha(range = alpha) +
                    guides(alpha = "none")
            } else if (label) {
                plot <- LabelClusters(
                    plot = plot,
                    id = ifelse(
                        test = is.null(x = cells.highlight),
                        yes = features[j],
                        no = 'highlight'
                    ),
                    geom = if (inherits(x = image.use, what = "STARmap")) {
                        'GeomPolygon'
                    } else if (use_geom_sf && plot_segmentations) {
                        'GeomSf'
                    } else if (use_geom_sf && !plot_segmentations) {
                        'GeomPoint'
                    } else {
                        'GeomSpatial'
                    },
                    repel = repel,
                    size = label.size,
                    color = label.color,
                    box = label.box,
                    position = "nearest"
                )
            }
            if (j == 1 && length(x = images) > 1 && !facet.highlight) {
                plot <- plot +
                    ggtitle(label = images[[image.idx]]) +
                    theme(plot.title = element_text(hjust = 0.5))
            }
            if (facet.highlight) {
                plot <- plot +
                    ggtitle(label = names(x = cells.highlight)[i]) +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    NoLegend()
            }
            if (use_geom_sf && plot_segmentations && !is.null(group.by)) {
                # Add legend guides to show filled squares next to labels when plotting segmentations
                plot <- plot + guides(fill = guide_legend(override.aes = list(alpha = 1, color = "black", linewidth = 0.2, size = 2)))
            }
            # Plot multiple images depending on keep.scale
            if (!(is.null(x = keep.scale)) && !inherits(x = data[, features[j]], "factor")) {
                plot <- suppressMessages(plot & scale_fill_gradientn(colors = SpatialColors(n = 100), limits = c(NA, max.feature.value)))
            }
            
            plots[[plot.idx]] <- plot
            plot.idx <- plot.idx + ncols
            if (cols.unset) {
                cols <- NULL
            }
        }
    }
    
    if (combine) {
        if (!is.null(x = ncol)) {
            return(patchwork::wrap_plots(plots = plots, ncol = ncol))
        }
        if (length(x = images) > 1) {
            return(patchwork::wrap_plots(plots = plots, ncol = length(x = images)))
        }
        return(patchwork::wrap_plots(plots = plots))
    }
    return(plots)
}
assignInNamespace("SpatialPlot",SpatialPlot1,ns="Seurat")
SpatialPlot = SpatialPlot1

SingleSpatialPlot1 <- function(
    data,
    image,
    cols = NULL,
    image.alpha = 1,
    image.scale = "lowres",
    pt.alpha = NULL,
    crop = TRUE,
    pt.size.factor = NULL,
    shape = 21,
    stroke = NA,
    col.by = NULL,
    alpha.by = NULL,
    cells.highlight = NULL,
    cols.highlight = c('#DE2D26', 'grey50'),
    geom = c('spatial', 'interactive', 'poly', 'sf'),
    na.value = 'grey50',
    plot_segmentations = FALSE
) {
  geom <- match.arg(arg = geom)
  if (!is.null(col.by) && !col.by %in% colnames(data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  
  col.by <- col.by %iff% paste0("`", col.by, "`")
  
  #Store unquoted col.by name for easier access
  col.by.clean <- gsub("`", "", col.by)
  
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size.factor,
      cols.highlight = cols.highlight[1],
      col.base = cols.highlight[2]
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), y = order))
    data <- data[order(data$ident), ]
  }
  plot <- ggplot(data = data, aes_string(
    x = colnames(data)[2],
    y = colnames(data)[1],
    fill = col.by,
    alpha = alpha.by
  ))
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      if (is.null(x = pt.alpha)) {
        plot <- plot + Seurat:::geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          image.scale = image.scale,
          crop = crop,
          shape = shape,
          stroke = stroke,
        )
      } else {
        plot <- plot + Seurat:::geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          image.scale = image.scale,
          crop = crop,
          shape = shape,
          stroke = stroke,
          alpha = pt.alpha
        )
      }
      plot + coord_fixed() + theme(aspect.ratio = 1)
    },
    'interactive' = {
      plot + Seurat:::geom_spatial_interactive(
        data = tibble::tibble(
          grob = list(
            GetImage(
              object = image,
              mode = 'grob'
            )
          )
        ),
        mapping = aes_string(grob = 'grob'),
        x = 0.5,
        y = 0.5
      ) +
        geom_point(mapping = aes_string(color = col.by)) +
        xlim(0, ncol(x = image)) +
        ylim(nrow(x = image), 0) +
        coord_cartesian(expand = FALSE)
    },
    'sf' = {
      
      # Validate image
      image.grob <- grid::rasterGrob(
        image@image,
        width = unit(1, "npc"),
        height = unit(1, "npc"),
        interpolate = FALSE
      )
      # Retrieve image dimensions for later use (flipping image)
      image.height <- dim(image@image)[1]
      image.width <- dim(image@image)[2]
      
      # Retrieve the sf data stored in the Visium V2 object 
      # Merge it with data dataframe which contains ident and gene expression information 
      sf.data = image@boundaries$segmentation@sf.data
      #Create sf object from data (POINTS), and extract xy
      data$cell <- rownames(data)
      data.sf <- sf::st_as_sf(data, coords = c("x", "y"), crs = NA)
      
      # Import pipe operator locally
      `%>%` <- magrittr::`%>%`
      
      data.coords <- data.sf %>%
        dplyr::mutate(x = sf::st_coordinates(.)[, 1],
               y = sf::st_coordinates(.)[, 2]) %>%
        sf::st_drop_geometry()
      
      #Merge with sf.data
      sf.merged <- sf.data %>%
        dplyr::left_join(data.coords, by = c("barcodes" = "cell"))
      sf.cleaned <- sf.merged %>% dplyr::filter(!is.na(x))
      
      #Extract centroids from VisiumV2, update sf.cleaned to match centroids
      if (!requireNamespace("sp", quietly = TRUE)) {
        stop("The 'sp' package is required but not installed.")
      }
      coordinates <- sp::coordinates
      coords <- coordinates(image@boundaries$centroids)
      barcodes <- image@boundaries$centroids@cells
      rownames(coords) <- barcodes
      
      # Find matching barcodes
      common_cells <- intersect(sf.cleaned$barcodes, rownames(coords))
      
      # Use match() to align indices
      match_idx <- match(common_cells, sf.cleaned$barcodes)
      coord_idx <- match(common_cells, rownames(coords))
      
      # Update x and y in sf.cleaned
      sf.cleaned$x[match_idx] <- coords[coord_idx, 1]
      sf.cleaned$y[match_idx] <- coords[coord_idx, 2]
      
      #Plot (currently independently of switch/case)
      if(!plot_segmentations){
        #If plot_segmentations FALSE, then plot just the centroids from sf data
        
        if (is.null(pt.alpha)) {
          #If pt.alpha not provided, then alpha parameter is derived from group/cluster data
          #Use alpha.by instead of pt.alpha
          geom_point_layer <- geom_point(
            shape = 21, 
            stroke = stroke, 
            size = pt.size.factor, 
            aes_string(fill = col.by, alpha = alpha.by)
          )
        } else {
          #If pt.alpha is indeed provided, then use that to define alpha
          geom_point_layer <- geom_point(
            shape = 21,
            stroke = stroke,
            size = pt.size.factor,
            aes_string(fill = col.by), 
            alpha = if (is.null(pt.alpha)) 1 else pt.alpha
          )
        }
        ggplot(sf.cleaned, aes(x = x, y = y)) +
          annotation_custom(
            grob = image.grob,
            xmin = 0,
            xmax = image.width,
            ymin = 0,
            ymax = image.height
          ) +
          geom_point_layer +
          coord_fixed() +
          xlab("x") +
          ylab("y") +
          theme_minimal()
        
      }else{
        #If plot_segmentations TRUE, then use geometry data stored in sf to plot cell polygons
        
        if (is.null(pt.alpha)) {
          #If pt.alpha not provided, then alpha parameter is derived from group/cluster data
          #Use alpha.by instead of pt.alpha
          geom_sf_layer <- geom_sf(
            data = sf.cleaned,
            aes_string(fill = col.by, alpha = alpha.by),
            color = "black",
            linewidth = stroke
          )
        } else {
          #If pt.alpha is indeed provided, then use that to define alpha
          geom_sf_layer <- geom_sf(
            data = sf.cleaned,
            aes_string(fill = col.by),
            alpha = if (is.null(pt.alpha)) 1 else pt.alpha,
            color = "black",
            linewidth = stroke
          ) 
        }
        ggplot() +
          annotation_custom(
            grob = image.grob,
            xmin = 0,
            xmax = image.width,
            ymin = 0,
            ymax = image.height
          ) +
          geom_sf_layer +
          scale_fill_viridis_d(option = "plasma") +
          coord_sf() +
          theme_void()
      }
    },
    'poly' = {
      data$cell <- rownames(x = data)
      data[, c('x', 'y')] <- NULL
      data <- merge(
        x = data,
        y = GetTissueCoordinates(
          object = image,
          qhulls = TRUE,
          scale = image.scale
        ),
        by = "cell"
      )
      plot + geom_polygon(
        data = data,
        mapping = aes_string(fill = col.by, group = 'cell')
      ) + coord_fixed() + theme_cowplot()
      
    },
    stop("Unknown geom, choose from 'spatial' or 'interactive'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = RColorBrewer::brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_fill_manual(values = colors, na.value = na.value)
    } else {
      data[[col.by.clean]] <- as.character(data[[col.by.clean]])
      vals <- unique(as.character(data[[col.by.clean]]))
      cols <- cols[names(cols) %in% vals]
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}
assignInNamespace("SingleSpatialPlot", SingleSpatialPlot1, ns="Seurat")
SingleSpatialPlot = SingleSpatialPlot1

CellsByImage <- function(object, images = NULL, unlist = FALSE) {
    images <- images %||% Images(object = object)
    cells <- sapply(
        X = images,
        FUN = function(x) {
            Cells(x = object[[x]])
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )
    if (unlist) {
        cells <- unname(obj = unlist(x = cells))
    }
    return(cells)
}

LabelClusters1 <- function(
        plot,
        id,
        clusters = NULL,
        labels = NULL,
        split.by = NULL,
        repel = TRUE,
        box = FALSE,
        geom = 'GeomPoint',
        position = "median",
        ...
) {
    xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom), use.names = TRUE)
    plot_data <- plot$data
    if (geom == "GeomSf") {
        # For sf, data is within the layers slot, not the data slot
        geom_layers <- which(sapply(plot$layers, function(layer) class(layer$geom)[1] == "GeomSf"))
        if (length(geom_layers) > 0 && !is.null(plot$layers[[geom_layers[1]]]$data)) {
            plot_data <- plot$layers[[geom_layers[1]]]$data
        }
    }
    if (!id %in% colnames(x = plot_data)) {
        stop("Cannot find variable ", id, " in plotting data")
    }
    if (!is.null(x = split.by) && !split.by %in% colnames(x = plot_data)) {
        warning("Cannot find splitting variable ", split.by, " in plotting data")
        split.by <- NULL
    }
    data <- plot_data[, c(xynames, id, split.by)]
    id_values <- if (inherits(data, "sf")) data[[id]] else data[, id]
    possible.clusters <- as.character(x = na.omit(object = unique(x = id_values)))
    groups <- clusters %||% possible.clusters
    if (any(!groups %in% possible.clusters)) {
        stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
    }
    pb <- ggplot_build(plot = plot)
    if (geom == 'GeomSpatial') {
        xrange.save <- layer_scales(plot = plot)$x$range$range
        yrange.save <- layer_scales(plot = plot)$y$range$range
        data[, xynames["y"]] = max(data[, xynames["y"]]) - data[, xynames["y"]] + min(data[, xynames["y"]])
        if (!pb$plot$plot_env$crop) {
            y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) - pb$layout$panel_params[[1]]$y.range
            data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
        }
    }
    data <- cbind(data, color = pb$data[[1]][[1]])
    labels.loc <- lapply(
        X = groups,
        FUN = function(group) {
            data.use <- if (inherits(data, "sf")) data[data[[id]] == group, , drop = FALSE] else data[data[, id] == group, , drop = FALSE]
            data.medians <- if (!is.null(x = split.by)) {
                do.call(
                    what = 'rbind',
                    args = lapply(
                        X = unique(x = data.use[, split.by]),
                        FUN = function(split) {
                            split_by_values <- if (inherits(data.use, "sf")) data.use[[split.by]] == split else data.use[, split.by] == split
                            split_data <- data.use[split_by_values == split, , drop = FALSE]
                            # Extract coordinates
                            if (inherits(split_data, "sf")) {
                                sf::st_agr(split_data) <- "constant" # Set attr-geom relationship to avoid warnings
                                coord_data <- data.frame(sf::st_coordinates(sf::st_centroid(split_data)))
                                names(coord_data) <- xynames[1:2]
                            } else {
                                coord_data <- data.use[data.use[, split.by] == split, xynames, drop = FALSE]
                            }
                            medians <- apply(
                                X = coord_data,
                                MARGIN = 2,
                                FUN = median,
                                na.rm = TRUE
                            )
                            medians <- as.data.frame(x = t(x = medians))
                            medians[, split.by] <- split
                            return(medians)
                        }
                    )
                )
            } else {
                # Extract coordinates
                if (inherits(data.use, "sf")) {
                    sf::st_agr(data.use) <- "constant"  # Set attr-geom relationship to avoid warnings
                    coord_data <- data.frame(sf::st_coordinates(sf::st_centroid(data.use)))
                    names(coord_data) <- xynames[1:2]
                } else {
                    coord_data <- data.use[, xynames, drop = FALSE]
                }
                as.data.frame(x = t(x = apply(
                    X = coord_data,
                    MARGIN = 2,
                    FUN = median,
                    na.rm = TRUE
                )))
            }
            data.medians[, id] <- group
            data.medians$color <- data.use$color[1]
            return(data.medians)
        }
    )
    if (position == "nearest") {
        labels.loc <- lapply(X = labels.loc, FUN = function(x) {
            # Handle sf data subsetting for nearest point calculation
            if (inherits(data, "sf")) {
                group.data <- data[as.character(data[[id]]) == as.character(x[3]), ]
                sf::st_agr(group.data) <- "constant"  # Set attr-geom relationship to avoid warnings
                group.data <- data.frame(sf::st_coordinates(sf::st_centroid(group.data)))
                names(group.data) <- xynames[1:2]
            } else {
                group.data <- data[as.character(x = data[, id]) == as.character(x[3]), ]
            }
            coord_matrix <- as.matrix(group.data[, 1:2])
            nearest.point <- RANN::nn2(data = coord_matrix, query = as.matrix(x = x[c(1,2)]), k = 1)$nn.idx
            x[1:2] <- coord_matrix[nearest.point, ]
            return(x)
        })
    }
    labels.loc <- do.call(what = 'rbind', args = labels.loc)
    # Safe handling of factor levels for sf data
    data_levels <- if (inherits(data, "sf")) levels(data[[id]]) else levels(data[, id])
    labels.loc[, id] <- factor(x = labels.loc[, id], levels = data_levels)
    labels <- labels %||% groups
    if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
        stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = unique(x = labels.loc[, id])), ").")
    }
    names(x = labels) <- groups
    for (group in groups) {
        labels.loc[labels.loc[, id] == group, id] <- labels[group]
    }
    if (box) {
        geom.use <- ifelse(test = repel, yes = ggrepel::geom_label_repel, no = geom_label)
        plot <- plot + geom.use(
            data = labels.loc,
            mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id, fill = id),
            show.legend = FALSE,
            ...
        )
    } else {
        geom.use <- ifelse(test = repel, yes = ggrepel::geom_text_repel, no = geom_text)
        plot <- plot + geom.use(
            data = labels.loc,
            mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
            show.legend = FALSE,
            ...
        )
    }
    # restore old axis ranges
    if (geom == 'GeomSpatial') {
        plot <- suppressMessages(expr = plot + coord_fixed(xlim = xrange.save, ylim = yrange.save))
    }
    return(plot)
}
assignInNamespace("LabelClusters", LabelClusters1, ns="Seurat")
LabelClusters = LabelClusters1

GetXYAesthetics1 <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
    geoms <- sapply(
        X = plot$layers,
        FUN = function(layer) {
            return(class(x = layer$geom)[1])
        }
    )
    # handle case where raster is set to True
    if (geom == "GeomPoint" && "GeomScattermore" %in% geoms){
        geom <- "GeomScattermore"
    }
    geoms <- which(x = geoms == geom)
    if (length(x = geoms) == 0) {
        stop("Cannot find a geom of class ", geom)
    }
    geoms <- min(geoms)
    if (plot.first) {
        # x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
        x <- rlang::as_label(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)
        # y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
        y <- rlang::as_label(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)
    } else {
        x <- rlang::as_label(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)
        y <- rlang::as_label(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)
    }
    # Handle GeomSf case where x/y are NULL because coordinates are in geometry
    if (geom == "GeomSf") {
        x <- "x"  # Default coordinate names for sf objects
        y <- "y"
    }
    return(list('x' = x, 'y' = y))
}
assignInNamespace("GetXYAesthetics", GetXYAesthetics1, ns="Seurat")
GetXYAesthetics = GetXYAesthetics1

SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))


################################################################################
