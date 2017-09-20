import errno    
import os

import statistics


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
#=pod

#=head2 writeCategory()

#Main function. Creates the HTML pages of the specified analysis.

#=cut

def writeCategory(outputdir, project_name, dataset_id, datatype, menu, d_samples):

#            Web::writeCategory($d_output, $d_htsnet, $d_rscript, $graphVizPath, $rh_datasets, $project, $study,
#                               "int", $rh_int_subnets_nwa, $rh_int_subnets_nwa, $rh_int_noa,  $rh_int_eda,
#                               $rh_int_go, $class, $sortByScore, $noHeatmap);
#    outputdir = shift;
#    itiPath = shift;
#    d_rscript = shift;
#    graphVizPath = shift;
#    rh_datasets = shift;
#    projName = shift;
#    dataset_id = shift;
#    network = shift;
#    rh_subnets_nwa = shift;
#    rh_subnetsOrMetas_nwa = shift;      # this variable equals $rh_subnets_nwa if we are not in the metasubnetwork case.
#    rh_noa = shift;
#    rh_eda = shift;
#    rh_GO = shift;
#    menu = shift;
#    sortByScore = shift;
#    noHeatmap = shift;
#    resHeading = shift;
#    noImage = shift;


#    # subnet part

#    my @a_snws = keys %$rh_subnets_nwa;
    level = 2

#    if(scalar(@a_snws) == 0){
#        subnetPageContent = "";
#        $subnetPageContent  = header($subnetPageContent, $projName, $dataset_id, $level, $network, "subnet", $menu, $resHeading);
#        $subnetPageContent .= "<h2>No subnetwork</h2>\n"
#        $subnetPageContent  = footer($subnetPageContent);

#        subnetPageName = "$dataset_id-$network-subnet.html";
#        open (SUBNET, ">$outputdir/$dataset_id/$network/$subnetPageName") || die "Cannot open $subnetPageName:\n$!";
#        print SUBNET $subnetPageContent;
#        close SUBNET;

#        genePageContent = "";
#        $genePageContent  = header($genePageContent, $projName, $dataset_id, $level, $network, "gene", $menu, $resHeading);
#        $genePageContent .= "<h2>No gene</h2>\n"
#        $genePageContent  = footer($genePageContent);

#        genePageName = "$dataset_id-$network-gene.html";
#        open (GENE, ">$outputdir/$dataset_id/$network/$genePageName") || die "Cannot open $genePageName:\n$!";
#        print GENE $genePageContent;
#        close GENE;

#        return;
#    }

#    my @a_datasets = sort keys %$rh_datasets;
#    
#    my @a_genes = keys %$rh_noa;

#    my @a_values;
#    for dataset(@a_datasets){
#        for gene(@a_genes){
#            val = $$rh_noa{$gene}{'Correlations'}{$dataset};
#            push(@a_values, abs($val)) if(defined $val);
#        }
#        $$rh_datasets{$dataset}{'Gradient_max'} = colorGradient(\@a_values);
#    }

#    if(defined $$rh_noa{$a_genes[0]}{'Highlights'}){
#        my %h_hlGenes;
#        for geneID(@a_genes){
#            for k(keys %{$$rh_noa{$geneID}{'Highlights'}}){
#                $h_hlGenes{$geneID} = 1 if(defined $$rh_noa{$geneID}{'Highlights'}{$k});
#            }
#        }
#        for snwID(keys %$rh_subnets_nwa){
#            nb = 0;
#            for geneID(keys %{$$rh_subnets_nwa{$snwID}{'Nodes'}}){
#                $nb++ if(defined $h_hlGenes{$geneID});
#            }
#            $$rh_subnets_nwa{$snwID}{'Highlighted_genes'} = $nb;
#        }

#        if($network eq 'itg'){
#            for snwID(keys %$rh_subnetsOrMetas_nwa){
#                nb = 0;
#                for geneID(keys %{$$rh_subnetsOrMetas_nwa{$snwID}{'Nodes'}}){
#                    $nb++ if(defined $h_hlGenes{$geneID});
#                }
#                $$rh_subnetsOrMetas_nwa{$snwID}{'Highlighted_genes'} = $nb;
#            }
#        }
#    }

#    # %h_dsSuperH('Transcriptome' => \%rh_datasets,
#    
#    my %h_dsSuperH;
#    $h_dsSuperH{'Transcriptome'} = $rh_datasets;
#    
#    # %h_dsLists('Transcriptome' => \@dataset_labels,
#    
#    my %h_dsLists;
#    $h_dsLists{'Transcriptome'} = \@a_datasets;
#    
    samplesdir = outputdir + "/" + dataset_id + "/" + datatype + "/samples";
    mkdir_p(samplesdir)
#    imageDir = "$outputdir/$dataset_id/$network/images";
#    mkdir $imageDir;

    samplePageContent = ""
    samplePageContent = sampleListPage(samplePageContent, project_name, dataset_id, datatype, menu, level, d_samples)

    samplePageName = dataset_id + "-" + datatype + "-sample.html";
    file = open(outputdir + "/" + dataset_id + "/" + datatype + "/" + samplePageName, "w") # || die "Cannot open $subnetPageName:\n$!";
    file.write(samplePageContent)
    file.close()


#    # gene part

#    geneDir = "$outputdir/$dataset_id/$network/genes";
#    mkdir $geneDir;

#    genePageContent = "";
#    genePageName = "$dataset_id-$network-gene.html";
#    $genePageContent = geneListPage($outputdir,      $genePageContent,       $projName, $dataset_id,$network, \%h_dsLists,
#                                    $rh_subnets_nwa, $rh_subnetsOrMetas_nwa, $rh_noa,    $menu,    $resHeading,
#                                    $rh_GO,          $level);

#    open (GENE, ">$outputdir/$dataset_id/$network/$genePageName") || die "Cannot open $genePageName:\n$!";
#    print GENE $genePageContent;
#    close GENE;

#    # GO part

#    if(defined $rh_GO){
#        goDir = "$outputdir/$dataset_id/$network/GO";
#        mkdir $goDir;

#        goPageContent = "";
#        goPageName = "$dataset_id-$network-go.html";
#        $goPageContent = goListPage($outputdir,      $goPageContent,        $projName, $dataset_id, $network, \%h_dsLists,
#                                    $rh_subnets_nwa, $rh_subnetsOrMetas_nwa, $rh_noa,    $menu,    $resHeading,
#                                    $rh_GO,          $level);

#        open (GO, ">$outputdir/$dataset_id/$network/$goPageName") || die "Cannot open $goPageName:\n$!";
#        print GO $goPageContent;
#        close GO;
#    }

#    #print $network." OK.\n"
#}

### SUBNETS ##

##=pod

##=head2 subnetListPage()

##Creates subnetwork pages, tables and images.

##=cut

def sampleListPage(filecontent, project_name, dataset_id, datatype, menu, level, d_samples):

#    outputdir = shift;
#    itiPath = shift;
#    d_rscript = shift;
#    graphVizPath = shift;
#    filecontent = shift;
#    rh_dsSuperH = shift;
#    rh_dsLists = shift;
#    projName = shift;
#    dataset_id = shift;
#    network = shift;
#    rh_subnets_nwa = shift;
#    rh_subnetsOrMetas_nwa = shift;
#    rh_noa = shift;
#    rh_eda = shift;
#    menu = shift;
#    resHeading = shift;
#    rh_GO = shift;
#    level = shift;
#    sortByScore = shift;
#    noImage = shift;
#    noHeatmap = shift;

    tab = "sample"

    filecontent = header("", project_name, dataset_id, datatype, menu, tab, level)

#    if($network eq "itg"){
#        metasubnetImage ($outputdir, $graphVizPath, $dataset_id, $network, $$rh_dsLists{'Transcriptome'}, $rh_subnets_nwa, $rh_subnetsOrMetas_nwa);

#        $filecontent = metasubnetBox($filecontent, $dataset_id, $network, $level);
#        $filecontent = metasubnetListBox($filecontent, $dataset_id, $network, $rh_subnetsOrMetas_nwa, $rh_noa, $level);
#    }
#    else {
    filecontent = sampleListBox(filecontent, dataset_id, datatype, level, d_samples)
#    }
    filecontent = footer(filecontent)

#    foreach subnetID (keys %$rh_subnetsOrMetas_nwa) {
#        subnetPage($outputdir,   $itiPath,   $d_rscript,    $graphVizPath, $projName, $dataset_id,      $network,
#                   $rh_dsSuperH, $rh_dsLists, $subnetID,     $rh_subnets_nwa, $rh_subnetsOrMetas_nwa,
#                   $rh_noa,      $rh_eda,     $menu,         $resHeading,     $rh_GO,
#                   $level, $sortByScore, $noImage, $noHeatmap);
#    }

    return filecontent


#=pod

#=head2 subnetListBox()

#Creates the table listing subnetworks and their information.

#=cut

def sampleListBox(filecontent, dataset_id, datatype, level, d_samples):
#    filecontent = shift;
#    dataset_id = shift;
#    network = shift;
#    rh_dsLists = shift;
#    rh_subnets_nwa = shift;
#    rh_subnetsOrMetas_nwa = shift;
#    rh_noa = shift;
#    level = shift;
#    goID = shift; ## defined if subnetlist from go page


    prefix = '../' * level

#    my @a_snws = keys %$rh_subnetsOrMetas_nwa;

    filecontent = filecontent + "<div class=\"box\">\n"
    filecontent = filecontent + "<h2>Samples (" + str(len(d_samples)) + ")</h2>\n"

    filecontent = filecontent + "<table class=\"sortable\" >\n"
    filecontent = filecontent + "<thead>\n"
    filecontent = filecontent + "<tr>\n"
    filecontent = filecontent + "<th><i class=\"fa fa-sort\"></i> Sample ID</th>\n"
    filecontent = filecontent + "<th><i class=\"fa fa-sort\"></i> Fastq file</th>\n"
    filecontent = filecontent + "<th><i class=\"fa fa-sort\"></i> Read number</th>\n"
    filecontent = filecontent + "<th><i class=\"fa fa-sort\"></i> Average read length</th>\n"
    filecontent = filecontent + "<th><i class=\"fa fa-sort\"></i> GC rate</th>\n"

#    for dataset(@{$$rh_dsLists{'Transcriptome'}}){
#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> $dataset score</th>\n"
#    }
#    
#    if(exists $$rh_subnets_nwa{$a_snws[0]}{'Merged_snw_best_score'}){
#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Merged snw best score</th>\n"
#    }

#    if(defined $$rh_subnets_nwa{$a_snws[0]}{'Seed_occurrence'}){
#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Seed occurrence</th>\n"
#    }
#    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Size</th>\n"

#    if(defined $$rh_subnets_nwa{$a_snws[0]}{'Highlighted_genes'}){
#        my @a_genes = keys %$rh_noa;
#        my @a_high  = keys %{$$rh_noa{$a_genes[0]}{'Highlights'}};
#        head    = "Highlighted genes";
#        if(scalar(@a_high) == 1){
#            $head = $a_high[0];
#            $head =~ s/_/ /g;
#        }

#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> $head</th>\n"
#    }

#    if(defined $goID){
#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Hypergeometric test</th>\n"
#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Bonferroni corrected p-value</th>\n"
#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Enrichment ratio</th>\n"
#        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Occurrences in subnetwork</th>\n"
#    }
    filecontent = filecontent + "</tr>\n"
    filecontent = filecontent + "</thead>\n"
    filecontent = filecontent + "<tbody>\n"

    for sample_id in d_samples.keys():

        (read_count, avg_read_length, gc_rate) = statistics.fastq_stats(d_samples[sample_id]["path"])

        filecontent = filecontent + "<tr>\n"
        filecontent = filecontent + "<td><a href=\"" + prefix + dataset_id + "/" + datatype + "/samples/" + dataset_id + "-" + datatype + "-sample-" + sample_id + ".html#sample\">" + sample_id + "</a></td>\n"
        filecontent = filecontent + "<td>" + d_samples[sample_id]["path"] + "</a></td>\n"
        filecontent = filecontent + "<td>" + str(read_count) + "</a></td>\n"
        filecontent = filecontent + "<td>" + str(avg_read_length) + "</a></td>\n"
        filecontent = filecontent + "<td>" + str(round(gc_rate, 1)) + "</a></td>\n"

#        score;
#        for dataset(@{$$rh_dsLists{'Transcriptome'}}){

#            $score = $$rh_subnets_nwa{$subnetID}{'Scores'}{$dataset};
#            if(defined $score){
#                $score = sprintf("%.3f", $score);
#            }
#            else{
#                $score = "-";
#            }

#            filecontent = filecontent + "<td>$score</td>\n"
#        }

#        if(defined $$rh_subnets_nwa{$subnetID}{'Merged_snw_best_score'}){
#            $score = $$rh_subnets_nwa{$subnetID}{'Merged_snw_best_score'};
#            if(defined $score){
#                $score = sprintf("%.3f", $score);
#            }
#            else{
#                $score = "-";
#            }
#            filecontent = filecontent + "<td>$score</td>\n"
#        }

#        if(defined $$rh_subnets_nwa{$subnetID}{'Seed_occurrence'}){
#            filecontent = filecontent + "<td>$$rh_subnets_nwa{$subnetID}{'Seed_occurrence'}</td>\n"
#        }
#        size = scalar(keys(%{$$rh_subnets_nwa{$subnetID}{'Nodes'}}));
#        filecontent = filecontent + "<td>$size</td>\n"

#        if(defined $$rh_subnets_nwa{$subnetID}{'Highlighted_genes'}){
#            filecontent = filecontent + "<td>$$rh_subnets_nwa{$subnetID}{'Highlighted_genes'}</td>\n"
#        }

#        if(defined $goID){
#            hyper = $$rh_subnetsOrMetas_nwa{$subnetID}{'GO_terms'}{$goID}{'hyper'};
#            corrpval = $$rh_subnetsOrMetas_nwa{$subnetID}{'GO_terms'}{$goID}{'corrected_pval'};
#            enrichRatio = $$rh_subnetsOrMetas_nwa{$subnetID}{'GO_terms'}{$goID}{'Enrichment_ratio'};
#            occ = $$rh_subnetsOrMetas_nwa{$subnetID}{'GO_terms'}{$goID}{'occurrences'};
#            filecontent = filecontent + "<td>$hyper</td>\n"
#            filecontent = filecontent + "<td>$corrpval</td>\n"
#            filecontent = filecontent + "<td>$enrichRatio</td>\n"
#            filecontent = filecontent + "<td>$occ</td>\n"
#        }
        filecontent = filecontent + "</tr>\n"

#    }
    filecontent = filecontent + "</tbody>\n"
    filecontent = filecontent + "</table>\n"

    filecontent = filecontent + "</div>\n"

    return filecontent

##=pod

##=head2 metasubnetListBox()

##Creates the table listing meta-subnetworks and their information.

##=cut

##sub metasubnetListBox {
##    filecontent = shift;
##    # rh_dsLists = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_meta_nwa = shift;
##    rh_noa = shift;
##    level = shift;

##    prefix = '../' x $level;

##    my @a_snws = keys %$rh_meta_nwa;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>".scalar(keys %{$rh_meta_nwa})." meta-subnetworks</h2>\n"

##    filecontent = filecontent + "<table class=\"sortable\" >\n"
##    filecontent = filecontent + "<thead>\n"
##    filecontent = filecontent + "<tr>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Subnetwork ID</th>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Size</th>\n"

##    if(defined $$rh_meta_nwa{$a_snws[0]}{'Highlighted_genes'}){
##        my @a_genes = keys %$rh_noa;
##        my @a_high  = keys %{$$rh_noa{$a_genes[0]}{'Highlights'}};
##        head    = "Highlighted genes";
##        if(scalar(@a_high) == 1){
##            $head = $a_high[0];
##            $head =~ s/_/ /g;
##        }

##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> $head</th>\n"
##    }

##    # for dataset(@{$$rh_dsLists{'Transcriptome'}}){
##    #     filecontent = filecontent + "<th><i class=\"icon-sort\"></i> $dataset score</th>\n"
##    # }

##    filecontent = filecontent + "</tr>\n"
##    filecontent = filecontent + "</thead>\n"
##    filecontent = filecontent + "<tbody>\n"
##    foreach subnetID (@a_snws) {
##        size = scalar(@{$$rh_meta_nwa{$subnetID}{'Interactors'}});
##        filecontent = filecontent + "<tr>
##<td><a href=\"$prefix$dataset_id/$network/subnets/$dataset_id-$network-subnet-$subnetID.html#subnet\">$subnetID</a></td>
##<td>$size</td>\n"

##        if(defined $$rh_meta_nwa{$subnetID}{'Highlighted_genes'}){
##            filecontent = filecontent + "<td>$$rh_meta_nwa{$subnetID}{'Highlighted_genes'}</td>\n"
##        }

##        # score;
##        # for dataset(@{$$rh_dsLists{'Transcriptome'}}){
##        #     $score = $$rh_meta_nwa{$subnetID}{'Scores'}{$dataset};
##        #     if ($network ne "itg" and defined $score){
##        #         $score = sprintf("%.3f", $score);
##        #     }else{
##        #         $score = "-";
##        #     }
##        #     filecontent = filecontent + "<td>$score</td>\n"
##        # }
##	
##        filecontent = filecontent + "</tr>\n"
##    }

##    filecontent = filecontent + "</tbody>\n"
##    filecontent = filecontent + "</table>\n"

##    filecontent = filecontent + "</div>\n"

##    return $filecontent;
##}

##=pod

##=head2 subnetPage()

##Creates the page of the specified subnetwork.

##=cut

#sub subnetPage {
##    outputdir = shift;
##    itiPath = shift;
##    d_rscript = shift;
##    graphVizPath = shift;
##    projName = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_dsSuperH = shift;
##    rh_dsLists = shift;
##    subnetID = shift;
##    rh_subnets_nwa = shift;
##    rh_subnetsOrMetas_nwa = shift;
##    rh_noa = shift;
##    rh_eda = shift;
##    menu = shift;
##    resHeading = shift;
##    rh_GO = shift;
##    level = shift;
##    sortByScore = shift;
##    noImage = shift;
##    noHeatmap = shift;

##    $level++;
#    tab = "subnet";

##    my @geneList = keys %{$$rh_subnetsOrMetas_nwa{$subnetID}{'Nodes'}};
##    my @edgeList = @{$$rh_subnetsOrMetas_nwa{$subnetID}{'Edges'}};

##    # heatmap image generation

##    unless(defined $noImage){
##        subnetImage($outputdir,                    $graphVizPath,          $dataset_id, $network, $$rh_dsSuperH{'Transcriptome'},
##                    $$rh_dsLists{'Transcriptome'}, $rh_subnetsOrMetas_nwa, $rh_noa,    $rh_eda,  $subnetID);

##    }

##    unless(defined $noHeatmap){
##        generateHeatmap($outputdir, $itiPath, $d_rscript, $dataset_id, $network, $$rh_dsSuperH{'Transcriptome'}, $rh_subnetsOrMetas_nwa, $rh_noa,  $subnetID, undef, $sortByScore);
##    }

##    # subnet page content

#    samplesdir = "$outputdir/$dataset_id/$network/subnets";

#    filecontent = "";
#    filename = "$dataset_id-$network-subnet-$subnetID.html";

#    $filecontent = header($filecontent, $projName, $dataset_id, $level,
#                          $network,     $tab,       $menu,
#                          $resHeading,  $rh_GO);

#    filecontent = filecontent + "<div class=\"title\">\n"
#    filecontent = filecontent + "<h1>".$subnetID."</h1>\n"
#    filecontent = filecontent + "</div>\n"

##    if($network eq "itg"){
##        $filecontent = metaSubBox($filecontent, $dataset_id, $network, $subnetID, $rh_dsLists,
##                                  $rh_subnets_nwa, $rh_subnetsOrMetas_nwa, $rh_noa, $level);
##    }
##    else {
##        $filecontent = pvalBox($filecontent, $subnetID, $rh_dsLists, $rh_subnetsOrMetas_nwa);
##    }

#    $filecontent = subnetBox($filecontent, $rh_dsSuperH,
#                             $dataset_id, $network, $subnetID,
#                             $level);
#    $filecontent = geneListBox($filecontent, $dataset_id,
#                               $network, $rh_dsLists, \@geneList,
#                               $rh_noa, $subnetID, $rh_subnets_nwa, $level);
##    $filecontent = interactionListBox($filecontent, $dataset_id,
##                                      $network, \@edgeList, $rh_noa,
##                                      $rh_eda, $level);
##    if(defined $rh_GO){
##        my @goList = sort{$$rh_subnetsOrMetas_nwa{$subnetID}{'GO_terms'}{$a}{'hyper'} <=> $$rh_subnetsOrMetas_nwa{$subnetID}{'GO_terms'}{$b}{'hyper'}}
##                     keys(%{$$rh_subnetsOrMetas_nwa{$subnetID}{'GO_terms'}});
##        $filecontent = goListBox($filecontent, $dataset_id,             $network,   \@goList,
##                                 $rh_GO,       $rh_subnetsOrMetas_nwa, $subnetID,  $level,
##                                 "subnetpage");
##    }

#    $filecontent = footer($filecontent);

#    # page writing

#    open (FILE, ">$samplesdir/$filename") || die "Cannot open $filename:\n$!";
#    print FILE $filecontent;
#    close FILE;
#}

##=pod

##=head2 pvalBox()

##Creates the table containing statistics (score and p-values) of the specified subnetwork.

##=cut

##sub pvalBox {
##    filecontent = shift;
##    subnetID = shift;
##    rh_dsLists = shift;
##    rh_subnets_nwa = shift;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>Statistics</h2>\n"

##    filecontent = filecontent + "<table >\n"
##    filecontent = filecontent + "<thead>\n"
##    filecontent = filecontent + "<tr>\n"
##    filecontent = filecontent + "<th>Dataset</th>\n"
##    filecontent = filecontent + "<th>Score</th>\n"

##    dataset = $$rh_dsLists{'Transcriptome'}[0];
##    my @a_pValHeads;
##    @a_pValHeads = keys %{$$rh_subnets_nwa{$subnetID}{'p-values'}} if(defined $$rh_subnets_nwa{$subnetID}{'p-values'});
##    my @a_dsPValHeads = sort(grep(/-$dataset$/, @a_pValHeads));
##    for dist(@a_dsPValHeads){
##        $dist =~ s/-.*$//;
##        $dist =~ s/dist//;
##        filecontent = filecontent + "<th>p-value $dist</th>\n"
##    }

##    filecontent = filecontent + "</tr>\n"
##    filecontent = filecontent + "</thead>\n"
##    filecontent = filecontent + "<tbody>\n"

##    for dataset(@{$$rh_dsLists{'Transcriptome'}}){
##        filecontent = filecontent + "<tr>\n"
##        filecontent = filecontent + "<td>$dataset</td>\n"

##        score = $$rh_subnets_nwa{$subnetID}{'Scores'}{$dataset};
##        if(defined $score){
##            $score = sprintf("%.3f", $score);
##        }
##        else{
##            $score = "-";
##        }
##        filecontent = filecontent + "<td>$score</td>\n"

##        my @a_datasetPVals = sort(grep(/-$dataset$/, @a_pValHeads));
##        foreach p (@a_datasetPVals){
##            pval = $$rh_subnets_nwa{$subnetID}{'p-values'}{$p};
##            if(defined $pval){
##                $pval = sprintf("%.2e", $pval);
##            }
##            else{
##                $pval = "-";
##            }
##            filecontent = filecontent + "<td>$pval</td>\n"
##        }
##        filecontent = filecontent + "</tr>\n"
##    }
##    
##    
##    filecontent = filecontent + "</tbody>\n"
##    filecontent = filecontent + "</table>\n"
##    filecontent = filecontent + "</div>\n"

##    return($filecontent);
##}

##=pod

##=head2 metaSubBox()

##Creates the table containing statistics (score and p-values) of the specified meta-subnetwork.

##=cut

##sub metaSubBox {
##    filecontent = shift;
##    dataset_id = shift;
##    network = shift;
##    metaID = shift;
##    rh_dsLists = shift;
##    rh_subnets_nwa = shift;
##    rh_meta_nwa = shift;
##    rh_noa = shift;
##    level = shift;

##    prefix = '../' x $level;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>Interacting subnets</h2>\n"

##    filecontent = filecontent + "<table >\n"
##    filecontent = filecontent + "<thead>\n"
##    filecontent = filecontent + "<tr>\n"
##    filecontent = filecontent + "<th>Subnetwork</th>\n"
##    filecontent = filecontent + "<th>Dataset</th>\n"
##    filecontent = filecontent + "<th>Score</th>\n"

##    dataset = $$rh_dsLists{'Transcriptome'}[0];
##    my @a_snws = @{$$rh_meta_nwa{$metaID}{'Interactors'}};
##    my @a_pValHeads = keys %{$$rh_subnets_nwa{$a_snws[0]}{'p-values'}};
##    my @a_dsPValHeads = sort(grep(/-$dataset$/, @a_pValHeads));
##    for dist(@a_dsPValHeads){
##        $dist =~ s/-.*$//;
##        $dist =~ s/dist//;
##        filecontent = filecontent + "<th>p-value $dist</th>\n"
##    }

##    filecontent = filecontent + "<th>Size</th>\n"

##    if(defined $$rh_subnets_nwa{$a_snws[0]}{'Highlighted_genes'}){
##        my @a_genes = keys %$rh_noa;
##        my @a_high  = keys %{$$rh_noa{$a_genes[0]}{'Highlights'}};
##        head    = "Highlighted genes";
##        if(scalar(@a_high) == 1){
##            $head = $a_high[0];
##            $head =~ s/_/ /g;
##        }

##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> $head</th>\n"
##    }

##    filecontent = filecontent + "</tr>\n"
##    filecontent = filecontent + "</thead>\n"
##    filecontent = filecontent + "<tbody>\n"

##    foreach subnetID (@a_snws){

##        netlink;
##        if ($network eq 'itg') { $netlink = $$rh_subnets_nwa{$subnetID}{'Type'}; }
##        else { $netlink = $network; }

##        size = scalar(keys(%{$$rh_subnets_nwa{$subnetID}{'Nodes'}}));
##        annot;
##        $annot = $$rh_subnets_nwa{$subnetID}{'Highlighted_genes'}
##            if(defined $$rh_subnets_nwa{$subnetID}{'Highlighted_genes'});

##        my @a_pValHeads = keys %{$$rh_subnets_nwa{$subnetID}{'p-values'}};
##        for dataset(@{$$rh_dsLists{'Transcriptome'}}){
##            score = $$rh_subnets_nwa{$subnetID}{'Scores'}{$dataset};
##            if(defined $score){
##                $score = sprintf("%.3f", $score);
##            }
##            else{
##                $score = "-";
##            }

##            filecontent = filecontent + "<tr>\n"
##            filecontent = filecontent + "<td><a href=\"$prefix$dataset_id/$netlink/subnets/$dataset_id-$netlink-subnet-$subnetID.html#subnet\">$subnetID</a></td>\n"
##            filecontent = filecontent + "<td>$dataset</td>\n"
##            filecontent = filecontent + "<td>$score</td>\n"

##            my @a_datasetPVals = sort(grep(/-$dataset$/, @a_pValHeads));
##            foreach p (@a_datasetPVals){
##                pval = $$rh_subnets_nwa{$subnetID}{'p-values'}{$p};
##                if(defined $pval){
##                    $pval = sprintf("%.2e", $pval);
##                }
##                else{
##                    $pval = "-";
##                }
##                filecontent = filecontent + "<td>$pval</td>\n"
##            }

##            filecontent = filecontent + "<td>$size</td>\n"
##            filecontent = filecontent + "<td>$annot</td>\n" if(defined $annot);

##            filecontent = filecontent + "</tr>\n"
##        }
##    }

##    filecontent = filecontent + "</tbody>\n"
##    filecontent = filecontent + "</table>\n"
##    filecontent = filecontent + "</div>\n"

##    return($filecontent);
##}

##=pod

##=head2 subnetBox()

##Creates the images links of the specified subnetwork.

##=cut

#sub subnetBox {
#    filecontent = shift;
#    rh_dsSuperH = shift;
#    dataset_id = shift;
#    network = shift;
#    subnetID = shift;
#    level = shift;

#    prefix = '../' x $level;

#    for dataset(sort(keys %{$$rh_dsSuperH{'Transcriptome'}})){
#        filecontent = filecontent + "<div class=\"box\">\n"
#        filecontent = filecontent + "<img class=\"caption\" title=\"$subnetID subnetwork with $dataset dataset\" src=\"$prefix$dataset_id/$network/images/subnetwork-$network-$subnetID-$dataset.png\" alt=\"$dataset-$subnetID subnetwork\"/>";
#        unless(defined $$rh_dsSuperH{'Transcriptome'}{$dataset}{'noHeatmap'}){
#            filecontent = filecontent + "<img class=\"caption\" title=\"$subnetID heatmap from $dataset dataset\" src=\"$prefix$dataset_id/$network/images/heatmap-$network-$subnetID-$dataset.png\" alt=\"$dataset-$subnetID heatmap\"/>";
#        }
#        filecontent = filecontent + "</div>\n"
#    }
#    return($filecontent);
#}

##=pod

##=head2 generateHeatmap()

##Creates heatmaps of the specified subnetwork genes.

##=cut

##sub generateHeatmap {
##    outputdir = shift;
##    itiPath = shift;
##    d_rscript = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_datasets = shift;
##    rh_subnetsOrMetas_nwa = shift;
##    rh_noa = shift;
##    subnetID = shift;
##    sortByScore = shift;

##    #print "Creating $subnetID heatmap(s)...\n"
##    imageDir = "$outputdir/$dataset_id/$network/images";

##    foreach dataset(keys %$rh_datasets) {
##        if(defined $$rh_datasets{$dataset}{'noHeatmap'}){
##            #print "No heatmap for dataset $dataset\n"
##            next;
##        }

##        expname = "heatmap-$network-$subnetID-$dataset.txt";
##        my ($pngname, $scKey, $dsKey);
##        
##        $scKey = 'Correlations'; $dsKey = $dataset; $pngname = "heatmap-$network-$subnetID-$dataset.png";

##        open (EXP, ">$imageDir/$expname") || die "Cannot open $expname:\n$!";

##        my @a_samples;

##        # In the case cond file is provided, samples are sorted by conditions
##        if(scalar(keys %{$$rh_datasets{$dataset}{'cond'}})) {
##            @a_samples = sort {$$rh_datasets{$dataset}{'cond'}{$a} cmp $$rh_datasets{$dataset}{'cond'}{$b}}
##                              keys %{$$rh_datasets{$dataset}{'cond'}};
##        }
##        else {
##            @a_samples = @{$$rh_datasets{$dataset}{'expe'}};
##        }

##        i = 0;
##        foreach expe (@a_samples){
##            ($i > 0) ? print EXP "\t$expe" : print EXP "$expe";
##            $i++;
##        }
##        print EXP "\n"

##        condNb = 0;
##        if(scalar(keys %{$$rh_datasets{$dataset}{'cond'}})) {
##            print EXP "Condition";
##            foreach expe (@a_samples){
##                print EXP "\t$$rh_datasets{$dataset}{'cond'}{$expe}";
##            }
##            print EXP "\n"
##            $condNb = 1;
##        }

##        elsif(scalar(keys %{$$rh_datasets{$dataset}{'mcond'}})) {
##            my @a_condNames = @{$$rh_datasets{$dataset}{'condNames'}};

##            for cond(@a_condNames){
##                print EXP $cond;

##                for expe(@a_samples){
##                    if(defined $$rh_datasets{$dataset}{'mcond'}{$expe}{$cond}){
##                        print EXP "\t$$rh_datasets{$dataset}{'mcond'}{$expe}{$cond}";
##                    }
##                    else{
##                        print EXP "\tUNDEF";
##                    }
##                }

##                print EXP "\n"
##            }

##            $condNb = scalar(@a_condNames);
##        }

##        my @a_geneIDs;

##        # Sort by gene score
##        if(defined $sortByScore){
##            my(%h_geneScores);
##            for g(keys(%{$$rh_subnetsOrMetas_nwa{$subnetID}{'Nodes'}})){
##                if(defined $$rh_noa{$g}{$scKey}{$dsKey}){
##                    $h_geneScores{$g} = $$rh_noa{$g}{$scKey}{$dsKey};
##                }
##                else{
##                    $h_geneScores{$g} = 0;
##                }
##            }
##            @a_geneIDs = sort{$h_geneScores{$a} <=> $h_geneScores{$b}} keys %h_geneScores;
##        }
##        # Sort by gene symbol
##        else{
##            @a_geneIDs = sort{$$rh_noa{$b}{'Gene_symbol'} cmp $$rh_noa{$a}{'Gene_symbol'}}
##                         keys(%{$$rh_subnetsOrMetas_nwa{$subnetID}{'Nodes'}});
##        }

##        foreach geneID (@a_geneIDs){
##            print EXP $$rh_noa{$geneID}{'Gene_symbol'};
##            foreach expe (@a_samples){

##                if(defined $$rh_datasets{$dataset}{'collapseddata'}{$geneID}{$expe}){
##                    print EXP "\t$$rh_datasets{$dataset}{'collapseddata'}{$geneID}{$expe}";
##                }
##                else{
##                    print EXP "\t0";    # annHeatmap2 does not support 'NAs'
##                }
##            }
##            print EXP "\n"
##        }

##        close(EXP);
##        system($d_rscript . "Rscript $itiPath/R/tableToHeatmap.R $imageDir/$expname $condNb $imageDir/$pngname");
##    }
##}

##=pod

##=head2 metasubnetBox()

##Creates the link of the meta-subnetworks connections image.

##=cut

##sub metasubnetBox {
##    filecontent = shift;
##    dataset_id = shift;
##    network = shift;
##    level = shift;

##    prefix = '../' x $level;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>Integrated network</h2>\n"
##    filecontent = filecontent + "<img class=\"caption\" title=\"...\" src=\"$prefix$dataset_id/$network/images/$network-metasubnet.png\" alt=\"metasubnets\"/>";
##    filecontent = filecontent + "</div>\n"

##    return($filecontent);
##}

##=pod

##=head2 subnetImage()

##Creates the image of the specified subnetwork.

##=cut

##sub subnetImage {
##    outputdir = shift;
##    graphVizPath = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_datasets = shift;
##    ra_datasets = shift;
##    rh_subnetsOrMetas_nwa = shift;
##    rh_noa = shift;
##    rh_eda = shift;
##    subnetID = shift;
##    

##    imageDir = "$outputdir/$dataset_id/$network/images";
##    dottype  = "digraph";

##    foreach dataset(@$ra_datasets) {

##        s_DOT = "$dottype G {\n node[style = filled];\n"
##        $s_DOT .= "overlap=false\n"
##        #$s_DOT .= "splines=true\n"
##        #node[style = filled];\n"

##        my($dotname, $pngname);
##	
##	$dotname = "subnetwork-$network-$subnetID-$dataset.dot";
##	$pngname = "subnetwork-$network-$subnetID-$dataset.png";

##        gradMax = $$rh_datasets{$dataset}{'Gradient_max'};

##        my @nodes = keys %{$$rh_subnetsOrMetas_nwa{$subnetID}{'Nodes'}};
##        my @a_edges = @{$$rh_subnetsOrMetas_nwa{$subnetID}{'Edges'}};

##        foreach geneID (@nodes){

##            shape = "\"triangle\"";
##            color = "\"black\"";
##            pen = 0.5;
##            periph = 1;

##            value;
##	    $value = $$rh_noa{$geneID}{'Correlations'}{$dataset};
##	    
##            geneSymbol = $$rh_noa{$geneID}{'Gene_symbol'};

##            if (defined $$rh_noa{$geneID}{'TF'}) { $shape = "\"box\""; }
##            else                                 { $shape = "\"egg\""; }

##            if($network eq "itg"){
##                if (Utils::in_array($geneID, $$rh_subnetsOrMetas_nwa{$subnetID}{'Connectors'})){
##                    $color = "\"black\"";
##                    $pen = 2;
##                }elsif (Utils::in_array($geneID, $$rh_subnetsOrMetas_nwa{$subnetID}{'Int_nodes'})){
##                    $color = "\"purple\"";
##                    $pen = 1;
##                }elsif (Utils::in_array($geneID, $$rh_subnetsOrMetas_nwa{$subnetID}{'Reg_nodes'})){
##                    $color = "\"blue\"";
##                    $pen = 1;
##                }
##            }
##            else{
##                if(defined $$rh_subnetsOrMetas_nwa{$subnetID}{'Merged_snw_seeds'}){
##                    $pen = 2 if(Utils::in_array($geneID, $$rh_subnetsOrMetas_nwa{$subnetID}{'Merged_snw_seeds'}));
##                }
##                elsif($$rh_subnetsOrMetas_nwa{$subnetID}{'Seed'} eq $geneID) {
##                    $pen = 2;
##                }
##            }

##            my @k;
##            @k = keys %{$$rh_noa{$geneID}{'Highlights'}} if(defined $$rh_noa{$geneID}{'Highlights'});
##            foreach k (@k){
##                if (defined $$rh_noa{$geneID}{'Highlights'}{$k}){
##                    $periph = 2;
##                }
##            }
##            if(defined $$rh_subnetsOrMetas_nwa{$subnetID}{'Recurrent_nodes'}){
##                if(Utils::in_array($geneID, $$rh_subnetsOrMetas_nwa{$subnetID}{'Recurrent_nodes'})){
##                    $periph = 2;
##                }
##                else{
##                    $periph = 1;
##                }
##            }

##            gradient;
##            if (defined $value and $value != 0 and $value =~ /^-?\d+\.?\d*$/ ) {
##                if(abs($value) >= $gradMax){
##                    $gradient = 255;
##                }
##                else{
##                    $gradient = sprintf('%.0f', abs($value)*(255/$gradMax));
##                }

##                if($value < 0) {
##                    $s_DOT .= "{node [shape=$shape, peripheries=$periph, color=$color, penwidth=$pen, fillcolor = \"#";
##                    $s_DOT .= sprintf("%02x",0);
##                    $s_DOT .= sprintf("%02x",150);
##                    $s_DOT .= sprintf("%02x",0);
##                    $s_DOT .= sprintf("%02x",$gradient);
##                    $s_DOT .= "\"] \"$geneSymbol\"}\n"
##                } elsif($value > 0) {
##                    $s_DOT .= "{node [shape=$shape, peripheries=$periph, color=$color, penwidth=$pen, fillcolor = \"#";
##                    $s_DOT .= sprintf("%02x",200);
##                    $s_DOT .= sprintf("%02x",0);
##                    $s_DOT .= sprintf("%02x",0);
##                    $s_DOT .= sprintf("%02x",$gradient);
##                    $s_DOT .= "\"] \"$geneSymbol\"}\n"
##                }
##            } else {
##                $s_DOT .= "{node [shape=$shape, peripheries=$periph, color=$color, penwidth=$pen, fillcolor = \"#";
##                $s_DOT .= sprintf("%02x",255);
##                $s_DOT .= sprintf("%02x",255);
##                $s_DOT .= sprintf("%02x",255);
##                $s_DOT .= "A0\"] \"$geneSymbol\"}\n"
##            }
##        }

##        foreach edgeID (@a_edges){
##            my @pair = split('-', $edgeID);
##            geneIDa = $pair[0];
##            geneIDb = $pair[1];
##            geneSymbolA = $$rh_noa{$geneIDa}{'Gene_symbol'};
##            geneSymbolB = $$rh_noa{$geneIDb}{'Gene_symbol'};

##            # next if (! defined      $$rh_eda{$edgeID}{'Direction'});
##            for it(keys %{$$rh_eda{$edgeID}}){
##                if( $$rh_eda{$edgeID}{$it}{'Direction'} eq "--" ){
##                    $s_DOT .= "subgraph \"$geneSymbolA$geneSymbolB\" {edge [dir=none, color=purple]"."\"$geneSymbolA\" -> \"$geneSymbolB\"}\n"     ## ! if meta
##                }
##                if( $$rh_eda{$edgeID}{$it}{'Direction'} eq "<>" || $$rh_eda{$edgeID}{$it}{'Direction'} eq ">" ){
##                    $s_DOT .= "edge [color=blue]\n"
##                    $s_DOT .= "\"$geneSymbolA\" -> \"$geneSymbolB\"\n"
##                }
##                if( $$rh_eda{$edgeID}{$it}{'Direction'} eq "<>" || $$rh_eda{$edgeID}{$it}{'Direction'} eq "<" ){
##                    $s_DOT .= "edge [color=blue]\n"
##                    $s_DOT .= "\"$geneSymbolB\" -> \"$geneSymbolA\"\n"
##                }
##            }
##        }
##        $s_DOT .= "}\n"

##        open (DOT, ">$imageDir/$dotname") || die "Cannot open $dotname:\n$!";
##        print DOT $s_DOT;
##        close DOT;

##        system("$graphVizPath/neato -Tpng $imageDir/$dotname -o $imageDir/$pngname");
##    }
##}

##=pod

##=head2 metasubnetImage()

##Creates the image of each meta-subnetwork.

##=cut

##sub metasubnetImage {

##    outputdir = shift;
##    graphVizPath = shift;
##    dataset_id = shift,
##    network = shift;
##    ra_datasets = shift;
##    rh_subnets_nwa = shift;
##    rh_meta_nwa = shift;

##    imageDir = "$outputdir/$dataset_id/$network/images";
##    dottype = "graph";
##    fillcolor = "\"cadetblue3\"";
##    color = "\"#00808064\"";
##    dotname = "$network-metasubnet.dot";
##    pngname = "$network-metasubnet.png";

##    s_DOT = "$dottype G {\nsplines=true\noverlap=false\nnode[style = filled];\n"

##    # print nodes
##    foreach metaID (keys %$rh_meta_nwa) {
##        foreach dataset(@$ra_datasets) {

##            shape = "\"octagon\"";
##            pen = 1;

##            $metaID =~ s/meta-//g;

##            type = $$rh_subnets_nwa{$metaID}{'Type'};
##            if ($type eq "int")     { $color = "\"purple\""; $fillcolor = "\"#9B30FF64\""; }
##            elsif($type eq "reg")   { $color = "\"blue\""; $fillcolor = "\"#0000ff64\""; }

##            $s_DOT .= "{node [shape=$shape, style=\"filled\", fillcolor=$fillcolor, color=$color, penwidth=$pen] \"$metaID\"}\n"
##        }
##    }

##    # print edges
##    foreach metaID (keys %$rh_meta_nwa) {
##        foreach dataset(@$ra_datasets) {
##            shape = "\"egg\"";
##            pen = 1;

##            my @nodes = @{$$rh_meta_nwa{$metaID}{'Interactors'}};

##            $metaID =~ s/meta-//g;

##            foreach node (@nodes) {

##                next if($node eq $metaID);

##                my @seed1 = split('-', $metaID);
##                seed1 = $seed1[-1];

##                my @seed2 = split('-', $node);
##                seed2 = $seed2[-1];

##                if ($seed1 < $seed2){
##                    my @connect = Metasubnet::overlapCount($metaID, $node, $rh_subnets_nwa, $rh_subnets_nwa);
##                    connect = scalar(@connect);
##                    for(i = 0; $i < $connect; $i++) {
##                        $s_DOT .= "\"$node\" -- \"$metaID\"\n"
##                    }
##                }
##            }
##        }
##    }

##    $s_DOT .= "}\n"

##    open (DOT, ">$imageDir/$dotname") || die "Cannot open $dotname:\n$!";
##    print DOT $s_DOT;
##    close DOT;

##    system("$graphVizPath/neato -Tpng $imageDir/$dotname -o $imageDir/$pngname");
##}

### GENES ##

##=pod

##=head2 geneListPage()

##Creates gene pages and tables.

##=cut

#sub geneListPage {
#    outputdir = shift;
#    filecontent = shift;
##    projName = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_dsLists = shift;
##    rh_subnets_nwa = shift;
##    rh_subnetsOrMetas_nwa = shift;
##    rh_noa= shift;
##    menu = shift;
##    resHeading = shift;
##    rh_GO = shift;
##    level = shift;

##    dataset = $$rh_dsLists{'Transcriptome'}[0];
##    my @geneList = sort{abs($$rh_noa{$b}{'Correlations'}{$dataset}) <=> abs($$rh_noa{$a}{'Correlations'}{$dataset})} keys %$rh_noa;

#    tab = "gene";

#    $filecontent = header($filecontent, $projName, $dataset_id, $level,
#                          $network,     $tab,       $menu,
#                          $resHeading,  $rh_GO);

#    $filecontent = geneListBox($filecontent, $dataset_id,
#                               $network, $rh_dsLists, \@geneList,
#                               $rh_noa, undef, undef, $level);

#    foreach geneID (@geneList) {
#        genePage($outputdir,      $projName, $dataset_id,           $network,    $geneID, $rh_noa,
#                 $rh_subnets_nwa, $rh_subnetsOrMetas_nwa, $rh_dsLists, $menu,   $resHeading,
#                 $rh_GO,          $level);
#    }

#    $filecontent = footer($filecontent);

#    return $filecontent;
#}

##=pod

##=head2 geneListBox()

##Creates the table listing specified genes and their information.

##=cut

#sub geneListBox {
#    filecontent = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_dsLists = shift;
##    ra_geneList = shift;
##    rh_noa = shift;
##    subnetID = shift;
##    rh_subnets_nwa = shift;
##    level = shift;

##    prefix = '../' x $level;

#    filecontent = filecontent + "<div class=\"box\">\n"
#    filecontent = filecontent + "<h2>Genes (".scalar(@$ra_geneList).")</h2>\n"

#    filecontent = filecontent + "<table class=\"sortable\" >\n"
#    filecontent = filecontent + "<thead>\n"
#    filecontent = filecontent + "<tr>\n"
#    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Gene Symbol</th>\n"
#    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Entrez Gene ID</th>\n"
#    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Frequency</th>\n"
##    seed;
##    comb = 0;
##    if(defined $rh_subnets_nwa){
##        $seed = $$rh_subnets_nwa{$subnetID}{'Seed'};
##        $comb = defined $$rh_subnets_nwa{$subnetID}{'Recurrent_nodes'};
##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Recurrent node</th>\n" if($comb);
##    }
##    for dataset(@{$$rh_dsLists{'Transcriptome'}}){
##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> $dataset gene score</th>\n"
##    }
##    
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Best subnetwork score</th>\n" if(defined $$rh_noa{$$ra_geneList[0]}{'Best_subnetwork_average_score'});

##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Degree</th>\n" if(defined $$rh_noa{$$ra_geneList[0]}{'Gene_degree'});

##    TF = exists($$rh_noa{$$ra_geneList[0]}{'TF'});
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Transcription factor</th>\n" if($TF);

##    my @high;
##    if (defined $$rh_noa{$$ra_geneList[0]}{'Highlights'}){
##        @high = keys %{$$rh_noa{$$ra_geneList[0]}{'Highlights'}};
##        # Warning! If we use @high the content is modified.
##        map {s/_/ /g; filecontent = filecontent + "<th><i class=\"icon-sort\"></i> $_</th>\n"}
##            keys %{$$rh_noa{$$ra_geneList[0]}{'Highlights'}};
##    }

#    filecontent = filecontent + "</tr>\n"
#    filecontent = filecontent + "</thead>\n"
#    filecontent = filecontent + "<tbody>\n"
##    for geneID (@$ra_geneList) {
##        freq = $$rh_noa{$geneID}{'Gene_frequency'};

##        sym;
##        if((defined $seed and $geneID eq $seed) ||
##           ($comb and Utils::in_array($geneID, $$rh_subnets_nwa{$subnetID}{'Merged_snw_seeds'}))){
##            $sym = "[ <a href=\"$prefix$dataset_id/$network/genes/$dataset_id-$network-gene-$geneID.html#gene\">$$rh_noa{$geneID}{'Gene_symbol'}</a> ]";
##        }
##        else{
##            $sym = "<a href=\"$prefix$dataset_id/$network/genes/$dataset_id-$network-gene-$geneID.html#gene\">$$rh_noa{$geneID}{'Gene_symbol'}</a>";
##        }

##        filecontent = filecontent + "<tr>\n"
##        $filecontent .=
##"<td>$sym</td>
##<td><a href=\"http://www.ncbi.nlm.nih.gov/gene/$geneID\" target=\"_blank\">$geneID</a></td>
##<td>$freq</td>";

##        if($comb){
##            if(Utils::in_array($geneID, $$rh_subnets_nwa{$subnetID}{'Recurrent_nodes'})){
##                filecontent = filecontent + "<td>Yes</td>";
##            } else{
##                filecontent = filecontent + "<td>-</td>";
##            }
##        }

##        for dataset(@{$$rh_dsLists{'Transcriptome'}}){
##            correlation = $$rh_noa{$geneID}{'Correlations'}{$dataset};

##            if(defined $correlation){
##                $correlation = sprintf("%.3f", $correlation);
##            }
##            else{
##                $correlation = "-";
##            }
##            if(defined $subnetID and defined $$rh_noa{$geneID}{'Present_in'}{$dataset}{$subnetID}){
##                filecontent = filecontent + "<td><strong>$correlation</strong></td>";
##            }
##            else{
##                filecontent = filecontent + "<td>$correlation</td>";
##            }
##        }
##	
##        score;
##        if(defined $$rh_noa{$geneID}{'Best_subnetwork_average_score'}){
##            $score = $$rh_noa{$geneID}{'Best_subnetwork_average_score'};
##            $score = sprintf("%.3f", $score);
##        }

##        filecontent = filecontent + "<td>$score</td>" if(defined $score);

##        degree = $$rh_noa{$geneID}{'Gene_degree'};
##        filecontent = filecontent + "<td>$degree</td>" if(defined $degree);

##        if($TF){
##            if(defined $$rh_noa{$geneID}{'TF'}){
##                filecontent = filecontent + "<td>TF</td>";
##            }
##            else{
##                filecontent = filecontent + "<td>-</td>";
##            }
##        }

##        foreach h (@high){
##            annot = $$rh_noa{$geneID}{'Highlights'}{$h};
##            if (defined $annot){
##                $annot =~ s/; /;<br>/g;
##                filecontent = filecontent + "<td>$annot</td>";
##            }
##            else{
##                filecontent = filecontent + "<td>-</td>";
##            }
##        }

##        filecontent = filecontent + "</tr>\n"
##    }
#    filecontent = filecontent + "</tbody>\n"
#    filecontent = filecontent + "</table>\n"

#    filecontent = filecontent + "</div>\n"

#    return $filecontent;
#}

##=pod

##=head2 interactionListBox()

##Creates the table listing specified interactions and their information.

##=cut

##sub interactionListBox {
##    filecontent = shift;
##    dataset_id = shift;
##    network = shift;
##    ra_interactionList = shift;
##    rh_noa = shift;
##    rh_eda = shift;
##    level = shift;

##    prefix = '../' x $level;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>Interactions (".scalar @$ra_interactionList.")</h2>\n"

##    filecontent = filecontent + "<table class=\"sortable\" >\n"
##    filecontent = filecontent + "<thead>\n"
##    filecontent = filecontent + "<tr>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Gene Symbol 1 </th>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Entrez Gene ID 1 </th>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Gene Symbol 2</th>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Entrez Gene ID 2</th>\n"

##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Type</th>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Direction</th>\n"
##    filecontent = filecontent + "<th>Origin databases / Sources</th>\n"
##    filecontent = filecontent + "</tr>\n"
##    filecontent = filecontent + "</thead>\n"
##    filecontent = filecontent + "<tbody>\n"

##    foreach edgeID (@{$ra_interactionList}) {
##        my @a_edge = split('-', $edgeID);
##        id1 = $a_edge[0];
##        id2 = $a_edge[1];

##        sym1 = $$rh_noa{$id1}{'Gene_symbol'};
##        sym2 = $$rh_noa{$id2}{'Gene_symbol'};

##        for type(keys %{$$rh_eda{$edgeID}}){
##            direction = $$rh_eda{$edgeID}{$type}{'Direction'};
##            source = $$rh_eda{$edgeID}{$type}{'Sources'};
##            my @a_source = split('\|\|', $source);

##            my @a_srcStr;
##            for dbS(@a_source){
##                my @a_dbS = split('::', $dbS);
##                db = $a_dbS[0];
##                so = $a_dbS[1];
##                if(defined $so){
##                    $so =~ s/;/, /g;
##                }
##                else{
##                    $so = "no annot";
##                }
##                push(@a_srcStr, "<strong>$db</strong>: $so");
##            }
##            srcStr = join(';<br>', @a_srcStr);

##            filecontent = filecontent + "<tr>\n"
##            $filecontent .=
##"<td><a href=\"$prefix$dataset_id/$network/genes/$dataset_id-$network-gene-$id1.html#gene\">".$sym1."</a></td>
##<td><a href=\"http://www.ncbi.nlm.nih.gov/gene/$id1\" target=\"_blank\">$id1</a></td>
##<td><a href=\"$prefix$dataset_id/$network/genes/$dataset_id-$network-gene-$id2.html#gene\">".$sym2."</a></td>
##<td><a href=\"http://www.ncbi.nlm.nih.gov/gene/$id2\" target=\"_blank\">$id2</a></td>
##<td>$type</td>
##<td>$direction</td>
##<td>$srcStr</td>\n"

##            filecontent = filecontent + "</tr>\n"
##        }
##    }
##    filecontent = filecontent + "</tbody>\n"
##    filecontent = filecontent + "</table>\n"

##    filecontent = filecontent + "</div>\n"

##    return $filecontent;
##}

##=pod

##=head2 genePage()

##Creates the page of the specified gene.

##=cut

#sub genePage {
#    outputdir = shift;
##    projName = shift;
##    dataset_id = shift;
##    network = shift;
##    geneID = shift;
##    rh_noa = shift;
##    rh_subnets_nwa = shift;
##    rh_subnetsOrMetas_nwa = shift;
##    rh_dsLists = shift;
##    menu = shift;
##    resHeading = shift;
##    rh_GO = shift;
##    level = shift;

##    $level++;
#    tab = "gene";
##    my @geneList = ($geneID);

##    geneDir = "$outputdir/$dataset_id/$network/genes";

##    my %h_geneSubnets = ();
##    foreach sub (@{$$rh_noa{$geneID}{'Subnets'}}) {
##        $h_geneSubnets{$sub} = $$rh_subnetsOrMetas_nwa{$sub};
##    }

#    filecontent = "";
#    filename = "$dataset_id-$network-$tab-$geneID.html";

#    $filecontent = header($filecontent, $projName, $dataset_id,
#                          $level, $network, $tab,
#                          $menu, $resHeading, $rh_GO);

#    filecontent = filecontent + "<div class=\"title\">\n"
#    filecontent = filecontent + "<h1>".$$rh_noa{$geneID}{'Gene_symbol'}."</h1>\n"
#    filecontent = filecontent + "</div>\n"

##    $filecontent = geneListBox($filecontent, $dataset_id,
##                               $network, $rh_dsLists, \@geneList,
##                               $rh_noa, undef, undef, $level);
##    $filecontent = geneBoxNCBI($filecontent, $geneID, $rh_noa);

##    if($network eq "itg"){
##        $filecontent = metasubnetListBox($filecontent, $dataset_id, $network, \%h_geneSubnets, $rh_noa, $level);
##    }
##    else{
##        $filecontent = subnetListBox($filecontent,    $dataset_id,      $network, $rh_dsLists,
##                                     $rh_subnets_nwa, \%h_geneSubnets, $rh_noa,  $level,      undef);
##    }

##    if(defined $rh_GO && defined $$rh_noa{$geneID}{'GO_terms'}){
##        my @goList = @{$$rh_noa{$geneID}{'GO_terms'}};
##        $filecontent = goListBox($filecontent, $dataset_id,
##                                 $network, \@goList, $rh_GO,
##                                 undef, undef, $level,
##                                 "genepage");
##    }

#    $filecontent = footer($filecontent);

#    open (FILE, ">$geneDir/$filename") || die "Cannot open $filename:\n$!";
#    print FILE $filecontent;
#    close FILE;
#}

##=pod

##=head2 geneBoxNCBI()

##Creates the NCBI link of the specified gene.

##=cut

##sub geneBoxNCBI{
##    filecontent = shift;
##    geneID = shift;
##    rh_noa = shift;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>Gene info in the <a href=\"http://www.ncbi.nlm.nih.gov\" target=\"_blank\">NCBI database</a>: <a href=\"http://www.ncbi.nlm.nih.gov/gene/$geneID\" target=\"_blank\">$$rh_noa{$geneID}{'Gene_symbol'}</a></h2>\n"

##    filecontent = filecontent + "</div>\n"

##    return($filecontent);
##}

#### GO ##

##=pod

##=head2 goListPage()

##Creates Gene Ontology pages and tables.

##=cut

##sub goListPage {
##    outputdir = shift;
##    filecontent = shift;
##    projName = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_dsLists = shift;
##    rh_subnets_nwa = shift;
##    rh_subnetsOrMetas_nwa = shift;
##    rh_noa= shift;
##    menu = shift;
##    resHeading = shift;
##    rh_GO = shift;
##    level = shift;

##    my @goList = keys %{$rh_GO};

##    tab = "go";
##    $filecontent = header($filecontent, $projName, $dataset_id, $level,
##                          $network,     $tab,       $menu,
##                          $resHeading,  $rh_GO);

##    $filecontent = goListBox($filecontent, $dataset_id,
##                             $network,     \@goList,     $rh_GO,
##                             undef,        undef,        $level,
##                             "golistpage");

##    foreach goID (@goList) {
##        goPage($outputdir,             $projName, $dataset_id, $network, $rh_dsLists, $rh_subnets_nwa,
##               $rh_subnetsOrMetas_nwa, $goID,      $rh_noa,  $menu,       $resHeading,
##               $rh_GO,                 $level);
##    }

##    $filecontent = footer($filecontent);

##    return $filecontent;
##}

##=pod

##=head2 goPage()

##Creates the page of the specified GO term.

##=cut

##sub goPage {
##    outputdir = shift;
##    projName = shift;
##    dataset_id = shift;
##    network = shift;
##    rh_dsLists = shift;
##    rh_subnets_nwa = shift;
##    rh_subnetsOrMetas_nwa = shift;
##    goID = shift;
##    rh_noa = shift;
##    menu = shift;
##    resHeading = shift;
##    rh_GO = shift;
##    level = shift;


##    $level++;
##    tab = "go";

##    my @goList = ($goID);

##    my %h_goSubnets = ();
##    foreach sub (@{$$rh_GO{$goID}{'subnets'}}) {
##        $h_goSubnets{$sub} = $$rh_subnetsOrMetas_nwa{$sub};
##    }

##    my @geneList = @{$$rh_GO{$goID}{'genes'}};

##    goDir = "$outputdir/$dataset_id/$network/GO";

##    filecontent = "";
##    filename = "$dataset_id-$network-$tab-$goID.html";

##    $filecontent = header($filecontent, $projName, $dataset_id, $level,
##                          $network,     $tab,       $menu,
##                          $resHeading,  $rh_GO);

##    filecontent = filecontent + "<div class=\"title\">\n"
##    filecontent = filecontent + "<h1>".$goID."</h1>\n"
##    filecontent = filecontent + "<h1>".$$rh_GO{$goID}{'term'}."</h1>\n"
##    filecontent = filecontent + "</div>\n"

##    $filecontent = goListBox($filecontent, $dataset_id, $network, \@goList,
##                             $rh_GO,     undef,        undef,      $level,   "gopage");
##    $filecontent = AmiGOBox($filecontent, $goID);
##    $filecontent = geneListBox($filecontent, $dataset_id,
##                               $network,     $rh_dsLists,  \@geneList,
##                               $rh_noa, undef, undef,      $level);
##    $filecontent = subnetListBox($filecontent, $dataset_id,
##                                 $network,     $rh_dsLists,  $rh_subnets_nwa,
##                                 \%h_goSubnets, $rh_noa, $level,       $goID);

##    $filecontent = footer($filecontent);

##    open (FILE, ">$goDir/$filename") || die "Cannot open $filename:\n$!";
##    print FILE $filecontent;
##    close FILE;
##}

##=pod

##=head2 AmiGOBox()

##Creates the AmiGO link of the specified GO term.

##=cut

##sub AmiGOBox{
##    filecontent = shift;
##    goID = shift;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>Term information on <a href=\"http://amigo.geneontology.org\" target=\"_blank\">AmiGO 2</a>: <a href=\"http://amigo.geneontology.org/amigo/term/$goID\" target=\"_blank\">$goID</a></h2>\n"
##    filecontent = filecontent + "</div>\n"

##    return($filecontent);
##}

##=pod

##=head2 goListBox()

##Creates the table listing specified GO terms and their information.

##=cut

##sub goListBox {
##    filecontent = shift;
##    dataset_id = shift;
##    network = shift;
##    ra_goList = shift;
##    rh_GO = shift;
##    rh_subnets_nwa = shift;
##    subnetID = shift;
##    level = shift;

##    prefix = '../' x $level;

##    filecontent = filecontent + "<div class=\"box\">\n"
##    filecontent = filecontent + "<h2>Related GO terms (".scalar(@{$ra_goList}).")</h2>\n"

##    filecontent = filecontent + "<p>Bold terms are significant according to the Benjamini-Hochberg correction (FDR = 5%).<br></p>" if(defined $subnetID);

##    filecontent = filecontent + "<table class=\"sortable\" >\n"
##    filecontent = filecontent + "<thead>\n"
##    filecontent = filecontent + "<tr>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Accession number</th>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Name</th>\n"
##    if(defined $subnetID){
##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Hypergeometric test</th>\n"
##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Bonferroni corrected p-value</th>\n"
##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Enrichment ratio</th>\n"
##        filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Occurrence in subnetwork</th>\n"
##    }
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Occurrences in all snw genes</th>\n"
##    filecontent = filecontent + "<th><i class=\"icon-sort\"></i> Occurrences in all int/reg genes</th>\n"

##    filecontent = filecontent + "</tr>\n"
##    filecontent = filecontent + "</thead>\n"
##    filecontent = filecontent + "<tbody>\n"
##    foreach goID (@{$ra_goList}) {

##        term = $$rh_GO{$goID}{'term'};
##        allSnwOcc = scalar @{$$rh_GO{$goID}{'genes'}};
##        totOcc = $$rh_GO{$goID}{'network_occurrences'};

##        filecontent = filecontent + "<tr>\n"
##        filecontent = filecontent + "<td><a href=\"$prefix$dataset_id/$network/GO/$dataset_id-$network-go-$goID.html#go\">".$goID."</a></td>";
##        if(defined $subnetID){
##            if(defined $$rh_subnets_nwa{$subnetID}{'GO_terms'}{$goID}{'BH_correction'} and
##               $$rh_subnets_nwa{$subnetID}{'GO_terms'}{$goID}{'BH_correction'} eq '*'){
##                $term = "<strong>$term</strong>";
##            }

##            hyper = $$rh_subnets_nwa{$subnetID}{'GO_terms'}{$goID}{'hyper'};
##            corrPVal = $$rh_subnets_nwa{$subnetID}{'GO_terms'}{$goID}{'corrected_pval'};
##            enrichRatio = $$rh_subnets_nwa{$subnetID}{'GO_terms'}{$goID}{'Enrichment_ratio'};
##            snwOcc = $$rh_subnets_nwa{$subnetID}{'GO_terms'}{$goID}{'occurrences'};
##            filecontent = filecontent + "<td>$term</td><td>$hyper</td><td>$corrPVal</td><td>$enrichRatio</td><td>$snwOcc</td>";
##        }
##        else{
##            filecontent = filecontent + "<td>$term</td>";
##        }
##        filecontent = filecontent + "<td>$allSnwOcc</td>";
##        filecontent = filecontent + "<td>$totOcc</td>";

##        filecontent = filecontent + "</tr>\n"
##    }
##    filecontent = filecontent + "</tbody>\n"
##    filecontent = filecontent + "</table>\n"

##    filecontent = filecontent + "</div>\n"

##    return $filecontent;
##}

# HOMEPAGE STUDY ##

#=pod

#=head2 homePage()

#Creates Home page of the specified analysis.

#=cut

def homePage(outputdir, project_name, dataset_id, datatype, description):

    level = 1
    tab = ""
    menu = "home"

    filecontent = header("", project_name, dataset_id, datatype, menu, tab, level)

    filecontent = studyBox(filecontent, dataset_id, description)

    filecontent = footer(filecontent)

    homePageName = dataset_id + "-home.html"

    mkdir_p(outputdir+ "/" + dataset_id)

    file = open (outputdir+ "/" + dataset_id + "/" + homePageName, "w") #|| die "Cannot open $homePageName:\n$!";
    file.write(filecontent)
    file.close()


#=pod

#=head2 studyBox()

#Writes the description of the specified analysis.

#=cut

def studyBox(filecontent, dataset_id, description):

    filecontent = filecontent + "<div class=\"box\">\n"
    filecontent = filecontent + "<h2>" + dataset_id + " analysis</h2>\n"
    filecontent = filecontent + "<p>" + description + "</p>\n"
    filecontent = filecontent + "</div>\n"

    return filecontent

### DOWNLOAD STUDY ##

##=pod

##=head2 downloadPage()

##Creates the Download page of the specified analysis.

##=cut

#sub downloadPage {
#    outputdir = shift;
#    dataset_id = shift;
#    menu = shift;
#    projName = shift;
#    resHeading = shift;
#    desc = shift;

#    level = 1;
#    tab = "";
#    network = "dl";

#    filecontent = header("", $projName, $dataset_id, $level,
#                             $network,     $tab,       $menu,
#                             $resHeading,  undef);

#    $filecontent = studyBox($filecontent, $dataset_id, $desc);
#    $filecontent = downloadBox($filecontent, $dataset_id, $menu,
#                               $level);

#    $filecontent = footer($filecontent);

#    downloadPageName = "$dataset_id-download.html";

#    open (DOWNLOAD, ">$outputdir/$dataset_id/$downloadPageName") || die "Cannot open $downloadPageName:\n$!";
#    print DOWNLOAD $filecontent;
#    close DOWNLOAD;
#}

def downloadPage(outputdir, project_name, dataset_id, datatype, description):

    level = 1
    tab = ""
    menu = "download"

    filecontent = header("", project_name, dataset_id, datatype, menu, tab, level)

    filecontent = downloadBox(filecontent, dataset_id, description)

    filecontent = footer(filecontent)

    downloadPageName = dataset_id + "-download.html"

    mkdir_p(outputdir+ "/" + dataset_id)

    file = open (outputdir+ "/" + dataset_id + "/" + downloadPageName, "w") #|| die "Cannot open downloadPageName:\n$!";
    file.write(filecontent)
    file.close()


#=pod

#=head2 downloadBox()

#Creates files links of the specified analysis.

#=cut

def downloadBox(filecontent, dataset_id, description):

    filecontent = filecontent + "<div class=\"box\">\n"
    filecontent = filecontent + "<h2>" + dataset_id + " analysis</h2>\n"
    filecontent = filecontent + "<p>" + description + "</p>\n"
    filecontent = filecontent + "</div>\n"

    return filecontent
#sub downloadBox {
#    filecontent = shift;
#    dataset_id = shift;
#    menu = shift;
#    level = shift;

#    prefix = '../' x $level;

#    filecontent = filecontent + "<div class=\"box\">";
#    filecontent = filecontent + "<h2>Data download</h2>";

#    if ($menu eq "itg"){
#        filecontent = filecontent + "<ul class=\"icons\">\n"
#        filecontent = filecontent + "<li><a href=\"$prefix$dataset_id/int/$dataset_id-int-subnet.html#subnet\"> Interactome</a></li>";
#        filecontent = filecontent + "</ul>\n"
#        filecontent = filecontent + "<table >\n"
#        filecontent = filecontent + "<tr><td>Download subnetworks file (*.nnf): </td><td><a href=\"$prefix$dataset_id/int/downloads/int-subnetworks.nnf\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download subnetwork statistics: </td><td><a href=\"$prefix$dataset_id/int/downloads/int-subnetworks_attributes-nwa.txt\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download gene list: </td><td><a href=\"$prefix$dataset_id/int/downloads/int-genes_attributes-noa.txt\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "</table><br/>\n"
#        filecontent = filecontent + "<ul class=\"icons\">\n"
#        filecontent = filecontent + "<li><a href=\"$prefix$dataset_id/reg/$dataset_id-reg-subnet.html#subnet\"> Regulome</a></li>";
#        filecontent = filecontent + "</ul>\n"
#        filecontent = filecontent + "<table >\n"
#        filecontent = filecontent + "<tr><td>Download subnetworks file (*.nnf): </td><td><a href=\"$prefix$dataset_id/reg/downloads/reg-subnetworks.nnf\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download subnetworks statistics: </td><td><a href=\"$prefix$dataset_id/reg/downloads/reg-subnetworks_attributes-nwa.txt\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download gene list: </td><td><a href=\"$prefix$dataset_id/reg/downloads/reg-genes_attributes-noa.txt\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "</table><br/>\n"
#        filecontent = filecontent + "<ul class=\"icons\">\n"
#        filecontent = filecontent + "<li><a href=\"$prefix$dataset_id/itg/$dataset_id-itg-subnet.html#subnet\"> Integrated Analysis - Meta Subnetworks</a></li>";
#        filecontent = filecontent + "</ul>\n"
#        filecontent = filecontent + "<table >\n"
#        filecontent = filecontent + "<tr><td>All results files (meta-subnetworks and their attributes): </td><td><a href=\"$prefix$dataset_id/itg/downloads/itg-results.zip\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download meta-subnetworks files (one *.nnf by meta-subnetwork): </td><td><a href=\"$prefix$dataset_id/itg/downloads/metasubnetworks.zip\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download meta-subnetworks information: </td><td><a href=\"$prefix$dataset_id/itg/downloads/metasubnets.txt\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download subnetworks list: </td><td><a href=\"$prefix$dataset_id/itg/downloads/subnets.txt\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Download genes information: </td><td><a href=\"$prefix$dataset_id/itg/downloads/meta_nodes.txt\" target=\"_blank\"><i class=\"icon-download-alt\"></i></a></td></tr>\n"
#        filecontent = filecontent + "</table>\n"
#    }
#    else{
#        filecontent = filecontent + "<table >\n"
#        filecontent = filecontent + "<tr><td>All results files (subnetworks and their attributes):</td><td><a href=\"$prefix$dataset_id/$menu/downloads/$menu-results.zip\" target=\"_blank\">Download All Files</a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Subnetworks files (one *.nnf by subnetwork):</td><td><a href=\"$prefix$dataset_id/$menu/downloads/$menu-subnetworks.zip\" target=\"_blank\">Download Nested Network Files</a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Subnetworks information:</td><td><a href=\"$prefix$dataset_id/$menu/downloads/$menu-subnetworks_attributes-nwa.txt\" target=\"_blank\">Download Subnetwork Attributes</a></td></tr>\n"
#        filecontent = filecontent + "<tr><td>Genes information:</td><td><a href=\"$prefix$dataset_id/$menu/downloads/$menu-genes_attributes-noa.txt\" target=\"_blank\">Download Node Attributes</a></td></tr>\n"
#        filecontent = filecontent + "</table><br/>\n"
#        filecontent = filecontent + "<ul class=\"icons\">\n"

#        filecontent = filecontent + "</ul>\n"
#    }

#    filecontent = filecontent + "</div>";

#    return $filecontent;
#}

## INDEX PROJECT ##

#=pod

#=head2 indexPage()

#Creates the Index page of the HTML report.

#=cut

def indexPage(outputdir, project_name, dataset_id):
#    outputdir = shift;
#    ra_studies = shift;
#    projName = shift;
#    resHeading = shift;

    level = 0
    datatype = ""
    tab = ""
    menu = ""
    filecontent = ""
    filename = "index.html"

    filecontent = header(filecontent, project_name, dataset_id, datatype, menu, tab, level)

    filecontent = indexBox(filecontent, project_name, dataset_id, level)

    filecontent = footer(filecontent)

    mkdir_p(outputdir)

    file = open (outputdir + "/" + filename, "w") # || die "Cannot open $filename:\n$!";
    file.write(filecontent)
    file.close()

#=pod

#=head2 indexBox()

#Creates analyses links in the Index page.

#=cut

def indexBox(filecontent, project_name, dataset_id, level):

    prefix = '../' * level;

    filecontent = filecontent + "<div class=\"box\">\n"
    filecontent = filecontent + "<h1>" + project_name + "</h1>\n"
    filecontent = filecontent + "<ul>\n"

#    foreach dataset_id (@$ra_studies){
    filecontent = filecontent + "<li><a href=\"" + prefix + dataset_id + "/" + dataset_id + "-home.html\">" + dataset_id + "</a></li>\n"
#    }

    filecontent = filecontent + "</ul>\n"
    filecontent = filecontent + "</div>\n"

    return filecontent

#=pod

#=head2 indexfooter()

#Creates the footer of the Index page.

#=cut

#sub indexfooter {
#    filecontent = shift;

#    filecontent = filecontent + "<div class=\"footer\"></div>\n"
#    filecontent = filecontent + "<div id=\"footer\">\n"
#    filecontent = filecontent + "<hr />\n"
#    filecontent = filecontent + "<img class=\"align-center\" src=\"KickStart/css/logo_cibi.png\" alt=\"cibi logo\"/><br/>\n"
#    filecontent = filecontent + "<br/>&copy; Copyright 2008-2015 <a href=\"http://cibi.marseille.inserm.fr/\">Cibi</a><br/>All Rights Reserved<br/>This website was built with <a href=\"http://www.99lime.com\">HTML KickStart</a>\n"

#    filecontent = filecontent + "</div>\n"
#    filecontent = filecontent + "</body>\n"
#    filecontent = filecontent + "</html>\n"

#    return $filecontent;
#}


## LAYOUT ##

#=pod

#=head2 header()

#Creates the standard header of web pages.

#=cut

def header(filecontent, project_name, dataset_id, datatype, menu, tab, level):

    prefix = '../' * level

    filecontent = filecontent + "<!doctype html>\n"
    filecontent = filecontent + "<html>\n"
    filecontent = filecontent + "<head>\n"
    filecontent = filecontent + "<title>" + project_name + "</title>\n"

#    filecontent = filecontent + "<link rel=\"icon\" type=\"image/png\" href=\"/home/rioualen/HTML-KickStart-master/css/sub-icon.png\" />\n"            ## !! abs path, to be parametrized

    filecontent = filecontent + "<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js\"></script>\n"
    filecontent = filecontent + "<script src=\"/home/rioualen/HTML-KickStart-master/js/kickstart.js\"></script>\n"                                      ## !! abs path, to be parametrized
    filecontent = filecontent + "<link rel=\"stylesheet\" href=\"/home/rioualen/HTML-KickStart-master/css/SnakeChunks.css\" media=\"all\" />\n"         ## !! abs path, to be parametrized

    filecontent = filecontent + "<meta http-equiv=\"Content-Type\" content=\"text/html;charset=utf-8\" />\n"

    filecontent = filecontent + "</head>\n"

    filecontent = filecontent + "<body>\n"

    if (level > 0):
        filecontent = menuBar(filecontent, dataset_id, datatype, menu, tab, level)

    filecontent = filecontent + "<div id=\"" + tab + "\" class=\"contenu\">\n"

    return filecontent


#=pod

#=head2 menuBar()

#Creates the menu of web pages.

#=cut

def menuBar(filecontent, dataset_id, datatype, menu, tab, level):

    prefix = '../' * level

    if (datatype == "chip"):
        chipclass = "\"current\""
        rnaclass = "\"\""
#        itgclass = "\"\""
        homeclass = "\"\""
        dlclass = "\"\""
#    }elsif ($network eq "int"){
#        $regclass = "\"\"";
#        $intclass = "\"current\"";
#        $itgclass = "\"\"";
#        $homeclass = "\"\"";
#        $dlclass = "\"\"";
    elif (datatype == "download"):
        chipclass = "\"\""
        rnaclass = "\"\""
        homeclass = "\"\""
        dlclass = "\"current\""
    elif (datatype == "home"):
        chipclass = "\"\""
        rnaclass = "\"\""
#        $itgclass = "\"\""
        homeclass = "\"current\""
        dlclass = "\"\""
#    }elsif ($network eq "itg"){
#        $regclass = "\"\"";
#        $intclass = "\"\"";
#        $itgclass = "\"current\"";
#        $homeclass = "\"\"";
#        $dlclass = "\"\"";
    
#    print(menu)

    filecontent = filecontent + "<div id=\"menu\">\n"
    filecontent = filecontent + "<h1><a href=\"" + prefix + "index.html\"><i class=\"fa fa-home\"></i></a> " + dataset_id + "</h1>        \n"

    filecontent = filecontent + "<ul class=\"menu\">\n"
    filecontent = filecontent + "<li class=" + homeclass + "><a href=\"" + prefix + dataset_id + "/" + dataset_id + "-home.html\"><i class=\"fa fa-home\"></i> Dataset</a></li>\n"
#    if($menu eq "itg"){
#        filecontent = filecontent + "<li class=$intclass><a href=\"$prefix$dataset_id/int/$dataset_id-int-subnet.html#subnet\">Interactome</a></li>\n"
#        filecontent = filecontent + "<li class=$regclass><a href=\"$prefix$dataset_id/reg/$dataset_id-reg-subnet.html#subnet\">Regulome</a></li>\n"
#        filecontent = filecontent + "<li class=$itgclass><a href=\"$prefix$dataset_id/itg/$dataset_id-itg-subnet.html#subnet\">Integrated analysis</a></li>\n"
#    }
#    else{
#    print(menu)
    filecontent = filecontent + "<li class=" + chipclass + "><a href=\"" + prefix + dataset_id + "/" + "chip" + "/" + dataset_id + "-" + "chip" + "-sample.html#sample\">ChIP-seq</a></li>\n"
#    }
    filecontent = filecontent + "<li class=" + dlclass + "><a href=\"" + prefix + dataset_id + "/" + dataset_id + "-download.html\"><i class=\"fa fa-download\"></i> Download</a></li>\n"
    filecontent = filecontent + "</ul>\n"

    if (level > 1):
        filecontent = tabsBar(filecontent, dataset_id, datatype, tab, prefix)

    filecontent = filecontent + "</div>\n"

    return filecontent


#=pod

#=head2 tabsBar()

#Creates tabs of web pages of the specified analysis.

#=cut

def tabsBar(filecontent, dataset_id, datatype, tab, prefix):

#    my ($subtab, $genetab, $gotab);
    if (tab == "peak"):
        sampletab = "\"" + prefix + dataset_id + + "/" + datatype + "/" + dataset_id + "-" + datatype + "-sample.html#sample\""
        peaktab = "\"#peak\""
    elif (tab == "sample"):
        sampletab = "\"#sample\"";
        peaktab = "\"" + prefix + dataset_id + "/" + datatype + "/" + dataset_id + "-" + datatype + "-peak.html#peak\""


    filecontent = filecontent + "<ul class=\"tabs left\">\n"
    filecontent = filecontent + "<li></li>\n"
    filecontent = filecontent + "<li><a href=" + sampletab + "> Samples</a></li>\n"
    filecontent = filecontent + "<li><a href=" + peaktab + "> Peaks</a></li>\n"
    filecontent = filecontent + "</ul>\n"

    return filecontent

#=pod

#=head2 footer()

#Creates the standard footer of web pages.

#=cut

def footer(filecontent):

#    filecontent = filecontent + "<div class=\"footer\"></div>\n"
    filecontent = filecontent + "<div id=\"footer\">\n"
    filecontent = filecontent + "<hr />This website was built with <a href=\"http://www.99lime.com\">HTML KickStart</a>\n"
    filecontent = filecontent + "</div>\n"
    filecontent = filecontent + "</div>\n"
    filecontent = filecontent + "</body>\n"
    filecontent = filecontent + "</html>\n"

    return filecontent


