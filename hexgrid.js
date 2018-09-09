
// const maxX = -35.2734, minX = -35.2848, maxY = 149.1360, minY = 149.1160;
const maxX = bounding_points[1][0], maxY = bounding_points[1][1], minX = bounding_points[0][0], minY = bounding_points[0][1]
const randomData = [];

for (i=0; i < 1000; i++) {
    const coordinates = [
        minX + Math.random() * (maxX - minX),
        minY + Math.random() * (maxY - minY)
    ];
    const properties = {
        a: 5 + Math.floor(Math.random() * 5),
        b: Math.floor(Math.random() * 5)
    };

    const marker = L.circleMarker(coordinates, {radius: 1, fillColor: 'black', fill: false, stroke: false});
    randomData.push( {
        marker: marker,
        properties: properties
    });
};

const rules = {
        cells: {
            "fillColor": {
                "method": "count",
                "attribute": "",
                "scale": "size",
                "range": ['#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'],
                "domain": [0,Math.min(numAgents/5.0,30)]
            },
            "color": "black",
            "fillOpacity": {
              "method": "count",
              "attribute": "",
              "scale": "size",
              "range": [0.0, 0.0174, 0.1570, 0.2346, 0.2792, 0.3277, 0.3800, 0.4363, 0.4963, 0.5603, 0.6, 0.6, 0.6, 0.6, 0.6], // [0, 0.1, 0.2, 0.4, 0.6]
              "domain": [0,Math.min(numAgents/2.0,30)]
            },
            "weight": 0
        },
        markers: {
            "color": "white",
            "weight": 2,
            "fillOpacity": 0.9,
            "fillColor": {
                "method": "count",
                "attribute": "",
                "scale": "continuous",
                "range": ["#ffffb2","#fecc5c","#fd8d3c","#e31a1c"]
            },
            "radius": {
                "method": "count",
                "attribute": "",
                "scale": "continuous",
                "range": [7, 17]
            }
        },
        texts: {
              "color": {
                "method": "count",
                "attribute": "",
                "scale": "quantile",
                "range": ["grey", "black"]
              },
              "fontSize": {
                "method": "count",
                "attribute": "",
                "scale": "continuous",
                "range": [12, 22]
              },
              "fontWeight": "bold",
              "text": {
                "method": "count",
                "attribute": ""
              },
              "anchorOffsetX": {
                "method": "count",
                "attribute": "",
                "scale": "continuous",
                "range": [3, 12]
              },
              "anchorOffsetY": {
                "method": "count",
                "attribute": "",
                "scale": "continuous",
                "range": [-5, -1]
              }
          }
    }

function creategrid() {
  const grid = L.regularGridCluster(
      {
          rules: rules,
          zoneSize: 7000*Math.pow(2,15-map.getZoom()), // at zoom 15, 10000 is about 312 m
          zoomHideGrid: 20,
          gridMode: 'hexagon',
          gridOrigin: L.latLng(minX-0.1, minY-0.1),
          gridEnd: L.latLng(maxX+0.1, maxY+0.1),
          showCells: true,
          showMarkers: false,
          showTexts: false
      });

  hexdata = []
  agentmap.agents.eachLayer(function(agent){
    const properties = {
        a: 5 + Math.floor(Math.random() * 5),
        b: Math.floor(Math.random() * 5)
    };
    hexdata.push({marker: agent, properties: properties})
  })

  grid.addLayers(hexdata)
  window.glayer = grid.addTo(map)
  return grid
}

grid = creategrid()
