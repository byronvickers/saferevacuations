
lineIntersect = require('@turf/line-intersect').default
lineDistance = require('@turf/point-to-line-distance').default
distance = require('@turf/distance').default
centroid = require('@turf/centroid').default
agentmaps = require('agentmaps').default
routing = require('agentmaps/src/routing')

////// dario code ////


bounding_points = [[-35.0626, 148.6739], [-35.5602, 149.5802]] // full act
let map = L.map("demo_map").fitBounds(bounding_points);
map_data.features.forEach((feature) => feature.properties.tags = { highway: true });

window.map = map;
L.tileLayer(
    "http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
    {
        attribution: "Thanks to <a href=\"http://openstreetmap.org\">OpenStreetMap</a> community",
    }
).addTo(map);

let street_options = {
  "color": "grey",
  "weight": 4,
  "opacity": .5
};

agentmap = L.A.agentmap(map);

let options_set = {
  standard: {
    "color": "blue",
    "weight": 1,
    "opacity": 0,
    "fillOpacity": 0,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 500,
    "depth": 500
  },
  suburb: {
    "color": "black",
    "weight": 1,
    "opacity": 0,
    "fillOpacity": 0,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 2000,
    "depth": 2000
  },
  person: {
    radius: .5,
    color: "green",
    fillColor: "green"
  }
};

let points_of_interest = [
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.237000,149.065000]], population: 97830, description: 'Belconnen'},
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.424400,149.088806]], population: 85968, description: 'Tuggeranong'},
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.168226,149.127139]], population: 72408, description: 'Gungahlin'},
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.280937,149.130005]], population: 54105, description: 'North Canberra'},
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.345200,149.095001]], population: 35377, description: 'Woden Valley'},
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.302684,149.122708]], population: 27663, description: 'South Canberra'},
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.341015,149.052521]], population: 23172, description: 'Weston Creek'},
  {type: 'suburb', unit_options: options_set.suburb, coords: [[-35.286015,149.063822]], population: 4759, description: 'Molonglo'}
];

agentmap.buildingify(bounding_points, map_data, street_options, options_set.standard, units_data);



/////

agentmap.units.eachLayer(function(unit){
  unit.options = Object.assign(unit.options, unit.feature.options)
  // unit.options.fillColor = unit.feature.options
  unit.setStyle()
})

////////

function destroygrid(grid) {
  if (map.hasLayer(grid)) {
    grid.unregister();
    map.removeLayer(grid);
  }
}

/* seeded js rng from https://stackoverflow.com/questions/521295/seeding-the-random-number-generator-in-javascript */
var seed = 1;
function random() {
    var x = Math.sin(seed++) * 10000;
    return x - Math.floor(x);
}

/* Durstenfeld shuffle from https://stackoverflow.com/questions/2450954/how-to-randomize-shuffle-a-javascript-array */
function shuffleArray(array) {
    for (let i = array.length - 1; i > 0; i--) {
        const j = Math.floor(random() * (i + 1));
        [array[i], array[j]] = [array[j], array[i]];
    }
}

function randomUnitAgentMaker(id){
	let index = this.agents.count();

  poitype = "suburb"

  unitset = this.units.getLayers().filter(unit => unit.feature.properties.poi_type === poitype)

  // let unit = this.units.getLayers()[Math.floor(agentmap.units.count() * Math.random())],
  let unit = unitset[Math.floor(unitset.length * Math.random())],
	unit_id = this.units.getLayerId(unit),
	center_point = centroid(unit.feature);
	center_point.properties.place = {"type": "unit", "id": unit_id},
	center_point.properties.layer_options = {radius: .5, color: "red", fillColor: "red"};

  return center_point;

}


window.numAgents = 50
idxarr = [...Array(numAgents).keys()];
shuffleArray(idxarr)
agentmap.agentify(numAgents, randomUnitAgentMaker);
// agentmap.agentify(numAgents, agentmap.seqUnitAgentMaker);

window.scatter = false
window.commandwaiting = false

congregate_idxs = [415, 301, 205]
agentmap.controller = function() {
    if ((agentmap.state.ticks === 0) || (window.commandwaiting)) {
        agentmap.agents.eachLayer(function(agent) {
              if(window.scatter){
                dest_index = Math.floor(agentmap.units.count() * Math.random())
                dest_unit = agentmap.units.getLayers()[dest_index]
              } else {
                closest_unit = null
                closest_dist = 9999999
                agentgeo = agent.toGeoJSON()
                congregate_idxs.forEach(function(idx){
                  unit = agentmap.units.getLayers()[idx]
                  unitgeo = centroid(unit.feature)
                  if(distance(agentgeo, unitgeo) < closest_dist){
                    closest_unit = unit
                    closest_dist = distance(agentgeo, unitgeo)
                  }
                })
                dest_unit = closest_unit
              }

              dest_unit_id = agentmap.units.getLayerId(dest_unit)
              dest_unit_center = dest_unit.getBounds().getCenter()

              if(agent.trip.path.length > 0){
                closest_unit = null
                closest_dist = 9999999
                agentgeo = agent.toGeoJSON()
                agentmap.units.eachLayer(function(unit){
                  unitgeo = centroid(unit.feature)
                  if(distance(agentgeo, unitgeo) < closest_dist){
                    closest_unit = unit
                    closest_dist = distance(agentgeo, unitgeo)
                  }
                })
                closest_unit_center = closest_unit.getBounds().getCenter()
                agent.setTravelToPlace(closest_unit_center, {type: "unit", id: agentmap.units.getLayerId(closest_unit)}, 25, true, true);
              }

              rndspeed = 4*(1 + random()/4)
              rnddelay = Math.floor(Math.random()*10000)
              // sendAgentWithDelay(agentmap.)
              try{
                agent.setTravelToPlace(dest_unit_center, {type: "unit", id: dest_unit_id}, rndspeed)
              }
              catch(error){
                // just leave them be, for now.
              }
        })
        window.commandwaiting = false
    }
    agentmap.agents.eachLayer(function(agent){
      agent.moveIt();
    })
    if (agentmap.state.ticks % 25 === 0){
      if(typeof grid !== 'undefined'){
        destroygrid(grid)
        grid = creategrid()
      }

    }
}
agentmap.run();
agentmap.setAnimationInterval(25)


// agentmap.streets.eachLayer(function(unit){unit.once("click", function(street){agentmap.streets.removeLayer(street)})})
// agentmap.streets.eachLayer(function(street){street.on("click", handleClick)})
map.on('click', function(e){
    clickCircle = new L.circle(e.latlng, 140*Math.pow(2, 15-map.getZoom()), {
      color: '#f07300',
      fillOpacity: 0.3,
      opacity: 0.5
    })
    .bindTooltip('🔥', { permanent: true, direction: "center", className: "labelstyle"})
    .addTo(map);

    // clickCircle.bringToFront()
    // window.c = clickCircle

    street_dists = []
    agentmap.streets.eachLayer(function(street){
      if(lineDistance(clickCircle.toGeoJSON(), street.feature)*1000 < clickCircle.getRadius()*0.8){ /* to get metres */
        agentmap.streets.removeLayer(street)
        street_id = agentmap.streets.getLayerId(street)

        agentmap.units.eachLayer(function(unit){
          if(unit.street_id === street_id){
            agentmap.units.removeLayer(unit)
          }
        });
      }
    })

    agentmap.streets.graph = routing.streetsToGraph(agentmap.streets),
  	agentmap.pathfinder = routing.getPathFinder(agentmap.streets.graph);

    window.commandwaiting = true
});

window.disperse = function(){
  window.scatter = true
  window.commandwaiting = true
}

window.congregate = function(){
  window.scatter = false
  window.commandwaiting = true
}

window.pausetoggle = function(){
  if(agentmap.state.paused){
    agentmap.run()
  } else{
    agentmap.pause()
  }
}
