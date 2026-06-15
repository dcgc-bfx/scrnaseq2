#' Add Plot Style
#'
#' Build a reusable set of ggplot2 theme and scale layers for consistent plotting
#' across figures.
#'
#' @param title Plot title. Default: NULL.
#' @param col Vector of colours to use for colour scales. Default: NULL.
#' @param fill Vector of colours to use for fill scales. Default: NULL.
#' @param legend_title Legend title to display or blank out. Default: NULL.
#' @param legend_position Legend position passed to `ggplot2::theme()`. Default: NULL.
#' @param xlab X-axis label to display or blank out. Default: NULL.
#' @param ylab Y-axis label to display or blank out. Default: NULL.
#' @param font_size Base font size for the theme. Default: 11.
#' @return A list of ggplot2 theme, label, and scale layers that can be added to a plot.
#' @note AI-assisted documentation
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

#' Generate Plot Captions
#'
#' Convert metric or feature names into readable plot captions, optionally splitting
#' paired names into comparison labels.
#'
#' @param plot_names Character vector or list of plot names to convert into captions.
#' @param assay_names Assay names to recognize and append in generated captions. Default: NULL.
#' @param split Separator used to split plot names into paired captions such as x vs y. Default: NULL.
#' @param capitalize Whether to capitalize the first letter of each caption. Default: TRUE.
#' @return A character vector of generated captions in the same order as `plot_names`.
#' @note AI-assisted documentation
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

#' Plot Barcode Quality Control
#'
#' Plot barcode-level quality-control metadata from a Seurat object, using violin
#' plots for numeric fields and stacked bar plots for categorical fields.
#'
#' @param sc Seurat v5 object containing barcode metadata.
#' @param qc Character vector of barcode metadata columns to plot.
#' @param filter Nested list of per-dataset filtering thresholds or allowed values for each QC field. Default: NULL.
#' @param assay Assay containing the barcodes to plot. Default: NULL, which uses the Seurat object's default assay.
#' @param log10 Whether to apply a log10 y-axis to numeric plots, or a regular expression selecting which metrics to transform. Default: FALSE.
#' @return A named list of ggplot2 objects, one plot per requested QC metric.
#' @note AI-assisted documentation
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

#' Plot Barcode QC Correlations
#'
#' Plot pairwise scatterplots for numeric barcode quality-control metrics and add
#' dataset-specific threshold lines where filters are provided.
#'
#' @param sc Seurat v5 object containing barcode metadata.
#' @param qc List of length-2 character vectors giving the barcode metadata column pairs to compare.
#' @param filter Nested list of per-dataset numeric filtering thresholds for QC fields. Default: NULL.
#' @param assay Assay containing the barcodes to plot. Default: NULL, which uses the Seurat object's default assay.
#' @return A named list of ggplot2 scatterplots, one plot per requested QC pair.
#' @note AI-assisted documentation
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

#' Plot Variable Features
#'
#' Create per-dataset scatterplots of variable-feature statistics and label the
#' top-ranked genes for the selected feature-selection method.
#'
#' @param sc Seurat v5 object containing variable-feature information.
#' @param method Variable-feature method to plot: `"vst"`, `"sct"`, or `"scran"`.
#' @param assay Assay to read variable-feature information from. Default: NULL, which uses the Seurat object's default assay.
#' @param top Number of top-ranked genes to label in each plot. Default: 10.
#' @return A named list of ggplot2 objects, one variable-feature plot per dataset.
#' @note AI-assisted documentation
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
                                       layer=paste("data", n, sep="."),
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

#' Plot Relative Log Expression
#'
#' Calculate relative log expression across sampled barcodes and summarize the
#' distribution for each dataset in a single QC plot.
#'
#' @param sc Seurat v5 object containing expression data.
#' @param assay Assay to use for the expression layer. Default: NULL, which uses the Seurat object's default assay.
#' @param layer Layer name or layer type to summarize, such as `"counts"` or `"data"`. Default: "counts".
#' @param nbarcodes Maximum number of barcodes to sample per dataset. Default: 500.
#' @param is_log Whether the selected layer is already log-transformed. Default: FALSE.
#' @return A ggplot2 object showing relative log expression summaries across sampled cells.
#' @note AI-assisted documentation
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

#' Create Spatial Dimensional Plots
#'
#' Dispatch to the appropriate Seurat spatial dimensional-plot function for each
#' requested image in a Seurat object.
#'
#' @param sc Seurat v5 object containing spatial images.
#' @param images Character vector of image names to plot. Default: NULL, which uses all images for the selected assay.
#' @param assay Assay used to discover images when `images` is NULL. Default: NULL, which uses the Seurat object's default assay.
#' @param ... Additional arguments passed to `Seurat::SpatialDimPlot()` or `Seurat::ImageDimPlot()`.
#' @return A named list of ggplot2 objects, one spatial dimensional plot per image.
#' @note AI-assisted documentation
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

#' Create Spatial Feature Plots
#'
#' Dispatch to the appropriate Seurat spatial feature-plot function for each
#' requested image in a Seurat object.
#'
#' @param sc Seurat v5 object containing spatial images.
#' @param images Character vector of image names to plot. Default: NULL, which uses all images for the selected assay.
#' @param assay Assay used to discover images when `images` is NULL. Default: NULL, which uses the Seurat object's default assay.
#' @param ... Additional arguments passed to `Seurat::SpatialFeaturePlot()` or `Seurat::ImageFeaturePlot()`.
#' @return A named list of ggplot2 objects, one spatial feature plot per image.
#' @note AI-assisted documentation
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

#' Build All Column Combinations Data Frame
#'
#' Reshape all pairwise column combinations from a matrix into a long-format data
#' frame suitable for downstream feature-grid style plotting.
#'
#' @param x Matrix with cells in rows and features or HTOs in columns.
#' @param cell_classification Vector of cell classifications matched to the matrix row names.
#' @return A data frame containing pairwise column combinations, values, and cell classifications for plotting.
#' @note AI-assisted documentation
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
