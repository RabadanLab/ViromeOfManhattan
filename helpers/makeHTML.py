import pandas as pd;

def generateHTML(reportLoc, namesLoc, outputLoc):
    """
    generateHTML: generate HTML report from report.taxon:
    
    reportLoc - location of report.taxon
    namesLoc - location of names.dmp file containing taxid to name
    outputLoc - location of file output folder
    """
    #reportLoc = "report.taxon.txt"
    #namesLoc = "/names.dmp"
    
    #outputLoc = "."
    #reportLoc = "~/report.taxon.txt"
    #namesLoc = "~/Downloads/taxdump/names.dmp"
    
    #load taxon report
    df = pd.read_table(reportLoc)
    df.applymap(str)
    df.columns = ['name','taxID','num_of_read','num_of_contigs','longest_contig','len_longest_contig','RPMH']
    sampleID = str(df['name'][0])
    
    ####
    #load canvasJS
    #    with open(canvasLoc,'r') as canvasjs:
    #        canvasString = canvasjs.read()
    #
    #    #load jquery
    #    with open(jqueryLoc,'r') as jqueryjs:
    #        jqueryString = jqueryjs.read()
    #
    #    #load DataTables
    #    with open(datatablesLoc,'r') as datatables:
    #        datatablesString = datatables.read()
    
    #load taxid to name
    nameDump = pd.read_table(namesLoc,sep="|",header=None)
    nameDump[1] = nameDump[1].str.strip()
    nameDump[2] = nameDump[2].str.strip()
    nameDump[3] = nameDump[3].str.strip()
    nameDump = nameDump.applymap(str)
    
    
    #generate canvasJS string for donut chart
    pd.options.mode.chained_assignment = None
    #donutString = ""
    df = df.applymap(str)
    for i in range(0,len(df.index)):
        df['taxID'][i] = df['taxID'][i].split(";")[0]
        tempname = nameDump[(nameDump[0] == df['taxID'][i]) & (nameDump[3] == 'scientific name')][[1]]
        if(not tempname.empty):
            df['name'][i] = tempname.to_string(header=False,index=False)
            #temp = "{y:"+str(df['RPMH'][i])+",label:'"+df['name'][i].replace("'","")+" -'},"
            #df['taxID'][i] = '<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id='+ str(df['taxID'][i]) +'">'+str(df['taxID'][i])+'</a>'
            #donutString=donutString+temp
    #donutString=donutString.rstrip(',')
    
    #workaround of column width to properly display taxid lnk
    #old_width = pd.get_option('display.max_colwidth')
    #pd.set_option('display.max_colwidth', -1)
    #tableString = df.to_html(index=False).replace("&lt;","<").replace("&gt;",">").replace('<table border="1" class="dataframe">','<table id="table1" border="1" class="dataframe">')
    #pd.set_option('display.max_colwidth', old_width)
    
    #generate output including JS packages as string in header
    htmloutput = '''
    <!DOCTYPE html>
    <html>
      <head>
        <meta charset="utf-8">
        <script src="https://d3js.org/d3.v4.min.js"></script>
      </head>
      <body>
        <div id="wrapper">
          <h1 style="margin-left:10%; font-size: 40px">Pandora Report</h1>
          <div id="chart"></div>
          <div id="table"></div>
        </div>
        <script>
          var chartDiv = document.getElementById("chart");
          var width = chartDiv.clientWidth -20;
          var height = chartDiv.clientHeight -20;
    
          var x = d3.scaleLinear().range([0,width]).domain([0,1000])
          var y = d3.scaleLinear().range([0,height]).domain([0,1000])
          var scale = d3.scaleLinear().range([0,width]).domain([0,1000])
    
          function donutChart() {
            var width,
                height,
                margin = {top: 10, right: 10, bottom: 10, left: 10},
                colour = d3.scaleOrdinal(d3.schemeCategory20c), // colour scheme
                variable, // value in data that will dictate proportions on chart
                category, // compare data by
                padAngle, // effectively dictates the gap between slices
                floatFormat = d3.format('.4r'),
                cornerRadius, // sets how rounded the corners are on each slice
                percentFormat = d3.format(',.2%');
    
            function chart(selection){
                selection.each(function(data) {
                    // generate chart
    
                    // ===========================================================================================
                    // Set up constructors for making donut. See https://github.com/d3/d3-shape/blob/master/README.md
                    var radius = Math.min(width, height) / 2;
    
                    // creates a new pie generator
                    var pie = d3.pie()
                        .value(function(d) { return floatFormat(d[variable]); })
                        .sort(null)
    
    
                    // contructs and arc generator. This will be used for the donut. The difference between outer and inner
                    // radius will dictate the thickness of the donut
                    var arc = d3.arc()
                        .outerRadius(radius * 0.8)
                        .innerRadius(radius * 0.5)
                        .cornerRadius(cornerRadius)
                        .padAngle(padAngle);
    
                    var bigArc = d3.arc()
                        .outerRadius(radius * 0.85)
                        .innerRadius(radius * 0.5)
                        .cornerRadius(cornerRadius)
                        .padAngle(padAngle);
    
                    // this arc is used for aligning the text labels
                    var outerArc = d3.arc()
                        .outerRadius(radius * 0.9)
                        .innerRadius(radius * 0.9);
                    // ===========================================================================================
    
                    // ===========================================================================================
                    // append the svg object to the selection
                    var svg = selection.append('svg')
                        .attr('width', width + margin.left + margin.right)
                        .attr('height', height + margin.top + margin.bottom)
                      .append('g')
                        .attr('transform', 'translate(' + width / 2 + ',' + height / 2 + ')');
                    // ===========================================================================================
    
                    // ===========================================================================================
                    // g elements to keep elements within svg modular
                    svg.append('g').attr('class', 'slices');
                    svg.append('g').attr('class', 'labelName');
                    svg.append('g').attr('class', 'lines');
                    // ===========================================================================================
    
                    // ===========================================================================================
                    // add and colour the donut slices
                    var path = svg.select('.slices')
                        .datum(data).selectAll('path')
                        .data(pie)
                      .enter().append('path')
                        .attr('fill', function(d) { return colour(d.data[category]); })
                        .attr('d', arc)
    
                    // ===========================================================================================
    
                    // ===========================================================================================
                    // add text labels
                    var label = svg.select('.labelName').selectAll('text')
                        .data(pie)
                      .enter().append('text')
                        .attr('dy', '.35em')
                        .html(function(d) {
                            // add "key: value" for given category. Number inside tspan is bolded in stylesheet.
                            return d.data[category] + ': <tspan>' + percentFormat(d.data[variable]) + '</tspan>';
                        })
                        .attr('transform', function(d) {
    
                            // effectively computes the centre of the slice.
                            // see https://github.com/d3/d3-shape/blob/master/README.md#arc_centroid
                            var pos = outerArc.centroid(d);
    
                            // changes the point to be on left or right depending on where label is.
                            pos[0] = radius * 0.95 * (midAngle(d) < Math.PI ? 1 : -1);
                            return 'translate(' + pos + ')';
                        })
                        .style('text-anchor', function(d) {
                            // if slice centre is on the left, anchor text to start, otherwise anchor to end
                            return (midAngle(d)) < Math.PI ? 'start' : 'end';
                        });
                    // ===========================================================================================
    
                    // ===========================================================================================
                    // add lines connecting labels to slice. A polyline creates straight lines connecting several points
                    var polyline = svg.select('.lines')
                        .selectAll('polyline')
                        .data(pie)
                      .enter().append('polyline')
                        .attr('points', function(d) {
    
                            // see label transform function for explanations of these three lines.
                            var pos = outerArc.centroid(d);
                            pos[0] = radius * 0.95 * (midAngle(d) < Math.PI ? 1 : -1);
                            return [arc.centroid(d), outerArc.centroid(d), pos]
                        });
                    // ===========================================================================================
    
                    // ===========================================================================================
                    // add tooltip to mouse events on slices and labels
                    d3.selectAll('.labelName text, .slices path').call(toolTip);
                    // ===========================================================================================
    
                    // ===========================================================================================
                    // Functions
    
                    // calculates the angle for the middle of a slice
                    function midAngle(d) { return d.startAngle + (d.endAngle - d.startAngle) / 2; }
    
                    // function that creates and adds the tool tip to a selected element
                    function toolTip(selection) {
    
                        // add tooltip (svg circle element) when mouse enters label or slice
                        selection.on('mouseenter', function (data) {
                            d3.select(this)
                              .transition()
                                .duration(200)
                                .attr('d',bigArc)
    
                            svg.append('text')
                                .attr('class', 'toolCircle')
                                .attr('dy', -15) // hard-coded. can adjust this to adjust text vertical alignment in tooltip
                                .html(toolTipHTML(data)) // add text to the circle.
                                .style('font-size', '.9em')
                                .style('text-anchor', 'middle'); // centres text in tooltip
    
                            svg.append('circle')
                                .attr('class', 'toolCircle')
                                .attr('r', radius * 0.45) // radius of tooltip circle
                                .style('fill', colour(data.data[category])) // colour based on category mouse is over
                                .style('fill-opacity', 0.35);
    
                        });
    
                        // remove the tooltip when mouse leaves the slice/label
                        selection.on('mouseout', function () {
                            d3.selectAll('.toolCircle').remove();
                            d3.select(this)
                              .transition()
                                .duration(100)
                                .attr('d',arc)
                        });
                    }
    
                    // function to create the HTML string for the tool tip. Loops through each key in data object
                    // and returns the html string key: value
                    function toolTipHTML(data) {
    
    
                        /*
                        var tip = '',
                            i   = 0;
    
                        for (var key in data.data) {
    
                            // if value is a number, format it as a percentage
                            var value = data.data[key];
    
                            // leave off 'dy' attr for first tspan so the 'dy' attr on text element works. The 'dy' attr on
                            // tspan effectively imitates a line break.
                            if (i === 0) tip += '<tspan x="0">' + key + ': ' + value + '</tspan>';
                            else tip += '<tspan x="0" dy="1.2em">' + key + ': ' + value + '</tspan>';
                            i++;
                        }
                        */
                        var tip = '<tspan x="0" dy="-0.5em" font-weight="bold">' + data.data.name + '</tspan>'
                          + '<tspan x="0" dy="1.2em">' + 'Percentage : ' + percentFormat(data.data[variable]) + '</tspan>'
                          + '<tspan x="0" dy="1.2em">' + 'Reads : ' + data.data.num_of_read + '</tspan>'
                          + '<tspan x="0" dy="1.2em">' + 'Contigs : ' + data.data.num_of_contigs + '</tspan>'
    
                        return tip;
                    }
                    // ===========================================================================================
    
                });
            }
    
            // getter and setter functions. See Mike Bostocks post "Towards Reusable Charts" for a tutorial on how this works.
            chart.width = function(value) {
                if (!arguments.length) return width;
                width = value;
                return chart;
            };
    
            chart.height = function(value) {
                if (!arguments.length) return height;
                height = value;
                return chart;
            };
    
            chart.margin = function(value) {
                if (!arguments.length) return margin;
                margin = value;
                return chart;
            };
    
            chart.radius = function(value) {
                if (!arguments.length) return radius;
                radius = value;
                return chart;
            };
    
            chart.padAngle = function(value) {
                if (!arguments.length) return padAngle;
                padAngle = value;
                return chart;
            };
    
            chart.cornerRadius = function(value) {
                if (!arguments.length) return cornerRadius;
                cornerRadius = value;
                return chart;
            };
    
            chart.colour = function(value) {
                if (!arguments.length) return colour;
                colour = value;
                return chart;
            };
    
            chart.variable = function(value) {
                if (!arguments.length) return variable;
                variable = value;
                return chart;
            };
    
            chart.category = function(value) {
                if (!arguments.length) return category;
                category = value;
                return chart;
            };
    
            return chart;
        }
    
        function barTable() {
          var width,
              cellheight,
              margin = {top: 10, right: 10, bottom: 10, left: 10},
              columns
    
          function chart(selection) {
            selection.each(function(data) {
              console.log(data)
              var svg = selection.append('svg')
                  .attr('width', '80%')
                  .attr('height', (data.length * cellheight) + margin.top + margin.bottom)
                  .attr('overflow','scroll')
                  .attr('transform', function(d,i){
                    return 'translate(' + scale(100) + ', 0)';
                  })
    
              var header = svg.append('g')
                .attr('class','header')
    
              //loop through object to generate columns
              for (var i in columns){
                var headerCol = header.append('g')
                  .attr('class','headerCol')
                  .attr('transform', 'translate(' + columns[i].x + ', 0)')
    
                headerCol.append('rect')
                  .attr('width',columns[i].width)
                  .attr('height',40)
                  .attr('rx',3)
                  .attr('ry',3)
                  .style('fill', 'rgb(253, 174, 107)')
                  .style('opacity',1)
    
                headerCol.append('text')
                  .attr('transform', function(d,i){
                    return 'translate(10, 30)';
                  })
                  .style('font-weight','bold')
                  .text(columns[i].name);
              }
    
              //draw table body
              var tablebody = svg.append('g')
                .attr('class','tablebody')
    
              var row = tablebody.selectAll('g')
                .data(data)
                .enter()
                .append('g')
                  .attr('class','row')
                  .attr('transform', function(d,i){
    
                    return 'translate(0, ' + ((i * cellheight) + 60)+ ')';
                  })
    
              var row_background = row.append('rect')
                .attr('transform', function(d,i){
                  return 'translate(' + scale(10) + ', -15)';
                })
                .attr('width','97%')
                .attr('height',20)
                .style('fill', 'rgb(253, 208, 162)')
                .style('opacity',0)
                .on('mouseover',function(d){
                  d3.select(this)
                    .style('opacity',0.9)
                })
                .on('mouseout',function(d){
                  d3.select(this)
                    .transition()
                    .duration(100)
                    .style('opacity',0)
                })
    
              var name_col = row.append('text')
                .attr('transform', function(d,i){
                  return 'translate(' + (10 + columns[0].x) + ', 0)';
                })
                .text(function(d,i){
                  return d.name
                });
    
              var taxid_col = row.append('a')
                .attr('transform', function(d,i){
                  return 'translate(' + (10 + columns[1].x) + ', 0)';
                })
                .attr('href', function(d){
                  return 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' + d.taxID
                })
                .attr('target','_blank')
                .append('text')
                  .text(function(d){
                    return d.taxID
                  });
    
              var longestcontigbox_col = row.append('rect')
                .attr('transform', function(d,i){
                  return 'translate(' + (5 + columns[2].x) + ', -15)';
                })
                .attr('height',20)
                .attr('width', function(d){
                  return d.contig_percentage * scale(90)
                })
                .attr('fill','rgb(198, 219, 239)')
                .attr('stroke','none')
    
              var longestcontig_col = row.append('text')
                .attr('transform', function(d,i){
                  return 'translate(' + (10 + columns[2].x) + ', 0)';
                })
                .text(function(d){
                  return d.len_longest_contig + ' bp'
                });
    
              var percentbox_col = row.append('rect')
                .attr('transform', function(d,i){
                  return 'translate(' + (5 + columns[3].x) + ', -15)';
                })
                .attr('height',20)
                .attr('width', function(d){
                  return d.percentage * scale(200)
                })
                .attr('fill','rgb(198, 219, 239)')
                .attr('stroke','none')
    
              var readnum_col = row.append('text')
                .attr('transform', function(d,i){
                  return 'translate(' + (10 + columns[3].x) + ', 0)';
                })
                .text(function(d){
                  return d.num_of_read
                });
    
              var rpmh_col = row.append('rect')
                .attr('transform', function(d,i){
                  return 'translate(' + (5 + columns[4].x) + ', -15)';
                })
                .attr('height',20)
                .attr('width', function(d){
                  return d.rpmh_percentage * scale(90)
                })
                .attr('fill','rgb(198, 219, 239)')
                .attr('stroke','none')
    
              var rpmh_col = row.append('text')
                .attr('transform', function(d,i){
                  return 'translate(' + (10 + columns[4].x) + ', 0)';
                })
                .text(function(d){
                  return d.RPMH
                });
    
            });
          }
          chart.width = function(value) {
              if (!arguments.length) return width;
              width = value;
              return chart;
          };
          chart.cellheight = function(value) {
              if (!arguments.length) return cellheight;
              cellheight = value;
              return chart;
          };
          chart.margin = function(value) {
              if (!arguments.length) return margin;
              margin = value;
              return chart;
          };
          chart.columns = function(value) {
              if (!arguments.length) return columns;
              columns = value;
              return chart;
          };
    
          return chart;
        }
    
        var donut = donutChart()
          .width(width)
          .height(scale(250))
          .cornerRadius(2) // sets how rounded the corners are on each slice
          .padAngle(0.015) // effectively dictates the gap between slices
          .variable('percentage')
          .category('name');
    
        var bartable = barTable()
          .width(width)
          .cellheight(22)
          .columns([
            {name: 'Name', width: scale(285), x:scale(10)},
            {name: 'TaxID', width: scale(95), x:scale(300)},
            {name: 'Longest contig', width: scale(95), x:scale(400)},
            {name: '# of reads', width: scale(195), x:scale(500)},
            {name: 'RPMH', width: scale(95), x:scale(700)},
          ])
    
        var data = d3.tsvParse(`''' + df.to_csv(sep="\t",index=False,doublequote=False) + '''`)
        //get total number of reads
        var total_num_reads = 0
        var longest_contig_length = 0
        var max_rpmh = 0
        for(var organism in data) {
          if(data[organism].num_of_read){
            total_num_reads += parseInt(data[organism].num_of_read)
          }
          if(parseInt(data[organism].len_longest_contig) > longest_contig_length){
            longest_contig_length = data[organism].len_longest_contig
          }
          if(parseInt(data[organism].RPMH) > max_rpmh){
            max_rpmh = data[organism].RPMH
          }
        }
    
        //filter out organisms over thresh
        for(var i in data){
          data[i].percentage = data[i].num_of_read / total_num_reads
          data[i].contig_percentage = data[i].len_longest_contig / longest_contig_length
          data[i].rpmh_percentage = data[i].RPMH / max_rpmh
        }
    
        data.sort(function(a, b) {
            return b.RPMH - a.RPMH ;
        });
    
        var filteredData = data.filter(function(d) {
          return d.percentage > 0.015
        })
        var underThresh = data.filter(function(d) {
          return d.percentage <= 0.015
        })
    
        //save organisms under thresh
        var readUnderThresh = 0
        for(var organism in underThresh) {
          if(underThresh[organism].num_of_read){
            readUnderThresh += parseInt(underThresh[organism].num_of_read)
          }
        }
    
        //add combined others to filteredData
        filteredData.push({
          name: "Others",
          num_of_read: readUnderThresh
        })
    
        for(var i in filteredData){
          filteredData[i].percentage = filteredData[i].num_of_read / total_num_reads
        }
    
        d3.select('#chart')
          .datum(filteredData) // bind data to the div
          .call(donut); // draw chart in div
    
        d3.select('#table')
          .datum(data)
          .call(bartable);
    
        </script>
      </body>
      <style>
        body {
          font-family: 'Roboto', sans-serif;
          color: #333333;
        }
        /*Styling for the lines connecting the labels to the slices*/
        polyline{
          opacity: .3;
          stroke: black;
          stroke-width: 1.5px;
          fill: none;
        }
    
        /* Make the percentage on the text labels bold*/
        .labelName tspan {
          font-style: normal;
          font-weight: 700;
        }
    
        /* In biology we generally italicise species names. */
        .labelName {
          font-size: 0.9em;
          font-style: italic;
        }
        #wrapper {
          display: table;
          height: 100vh;
        }
        #chart {
          display: table-row;
          left: 0px;
          right: 0px;
          top: 0px;
          bottom: 0px;
          z-index: 1;
        }
        #table {
          display:table-row;
          height:100%;
          overflow-y: scroll;
          font-size: 15px;
        }
      </style>
    </html>
    '''
    
    #write to html file
    with open(outputLoc+'/report.taxon.html', 'w') as f:
        f.write(htmloutput)
    
    
    
    
