function filterByGene(gene) {
	filtered = [];
	for (index in cancermine) {
		if (cancermine[index]['gene'] == gene)
			filtered.push(cancermine[index])
	}
	return gene;
}

//var data = [trace1, trace2];
//var layout = {barmode: 'stack'};

function test() {
	console.log('Hello');
}

selectedGene = 'EGFR1';

function update(selectedGene) {

	headerElement = document.getElementById('header');
	headerElement.innerHTML = 'Summary for gene: ' + selectedGene;

	driverLabels = []
	driverCounts = []
	oncogeneLabels = []
	oncogeneCounts = []
	suppressorLabels = []
	suppressorCounts = []
	
	totalDriverCounts = 0
	totalOncogeneCounts = 0
	totalSuppressorCounts = 0
	for (index in cancermine)
	{
		if (cancermine[index]['gene'] == selectedGene)
		{
			if (cancermine[index]['relationtype'] == 'Driver') {
				driverLabels.push(cancermine[index]['cancer']);
				driverCounts.push(cancermine[index]['count']);
				totalDriverCounts += cancermine[index]['count']
			} else if (cancermine[index]['relationtype'] == 'Oncogene') {
				oncogeneLabels.push(cancermine[index]['cancer']);
				oncogeneCounts.push(cancermine[index]['count']);
				totalOncogeneCounts += cancermine[index]['count']
			} else if (cancermine[index]['relationtype'] == 'Tumor Suppressor') {
				suppressorLabels.push(cancermine[index]['cancer']);
				suppressorCounts.push(cancermine[index]['count']);
				totalSuppressorCounts += cancermine[index]['count']
			}
		}
	}

	var data = [
		{x: driverLabels,y: driverCounts,name: 'Driver',type: 'bar'},
		{x: oncogeneLabels,y: oncogeneCounts,name: 'Oncogene',type: 'bar'},
		{x: suppressorLabels,y: suppressorCounts,name: 'Tumor Suppressor',type: 'bar'},
	];
	var layout = {
		height: 400,
		barmode: 'stack'
	};
	Plotly.newPlot('barchart', data, layout);
	
	var data = [{
		values: [totalDriverCounts, totalOncogeneCounts, totalSuppressorCounts],
		labels: ['Driver', 'Oncogene', 'Tumor Suppressor'],
		type: 'pie'
	}];
	
	var layout = {
	  height: 400,
	  width: 500,
	};
	
	Plotly.newPlot('piechart', data, layout);

	$(function () {
		$('#table').bootstrapTable({
			data: cancermine
		});
	});
}


update('EGFR1');

var input = document.getElementById("mysearch");
new Awesomplete(input, {
list: genes
});
document.querySelector('#mysearch').addEventListener('awesomplete-selectcomplete', function(evt){
  update(this.value);
})
