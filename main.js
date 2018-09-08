
lineIntersect = require('@turf/line-intersect').default
lineDistance = require('@turf/point-to-line-distance').default
distance = require('@turf/distance').default
centroid = require('@turf/centroid').default
agentmaps = require('agentmaps').default
routing = require('agentmaps/src/routing')

window.map = L.map("demo_map").fitBounds(bounding_points);
L.tileLayer(
    "http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
    {
        attribution: "Thanks to <a href=\"http://openstreetmap.org\">OpenStreetMap</a> community",
    }
).addTo(map);

map_data.features = map_data.features.filter(function(feature){return feature.geometry.type === "LineString" && feature.properties.tags.highway && (["trunk", "primary", "secondary", "tertiary", "turning_circle", "unclassified"].indexOf(feature.properties.tags.highway) >= 0) })
// map_data.features = map_data.features.filter(function(feature){return feature.geometry.type === "LineString" && feature.properties.tags.highway && (["trunk", "primary"].indexOf(feature.properties.tags.highway) >= 0) })

function destroygrid(grid) {
  if (map.hasLayer(grid)) {
    grid.unregister();
    map.removeLayer(grid);
  }
}

let street_options = {
 "color": "yellow",
 "weight": 4,
 "opacity": .5
};
let unit_options = {
 front_buffer: 6,
 side_buffer: 3,
 length: 70,
 depth: 18
}
agentmap = L.A.agentmap(map);
agentmap.buildingify(bounding_points, map_data, street_options, unit_options, units_data);

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
	// if (index > this.units.getLayers().length - 1) {
	// 	throw new Error("randomUnitAgentMaker cannot accommodate more agents than there are units.");
	// }

	// let unit = this.units.getLayers()[idxarr[index]],
  let unit = this.units.getLayers()[Math.floor(agentmap.units.count() * Math.random())],
	unit_id = this.units.getLayerId(unit),
	center_point = centroid(unit.feature);
	center_point.properties.place = {"type": "unit", "id": unit_id},
	center_point.properties.layer_options = {radius: .5, color: "red", fillColor: "red"};
  return center_point;

}

window.numAgents = 20
idxarr = [...Array(numAgents).keys()];
shuffleArray(idxarr)
agentmap.agentify(numAgents, randomUnitAgentMaker);
// agentmap.agentify(numAgents, agentmap.seqUnitAgentMaker);

window.scatter = false
window.commandwaiting = false

agentmap.controller = function() {
    if ((agentmap.state.ticks === 0) || (window.commandwaiting)) {
        agentmap.agents.eachLayer(function(agent) {
              if(window.scatter){
                dest_index = Math.floor(agentmap.units.count() * Math.random())
              } else {
                dest_index = 0
              }
              dest_unit = agentmap.units.getLayers()[dest_index]
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
                agent.setTravelToPlace(closest_unit_center, {type: "unit", id: agentmap.units.getLayerId(closest_unit)}, 1, true, true);
              }

              rndspeed = 1 + random()/4
              rnddelay = Math.floor(Math.random()*10000)
              // sendAgentWithDelay(agentmap.)
              agent.scheduleTrip(dest_unit_center, {type: "unit", id: dest_unit_id}, rndspeed)
        })
        window.commandwaiting = false
    }
    agentmap.agents.eachLayer(function(agent){
      agent.moveIt();
    })
    if (agentmap.state.ticks % 100 === 0){
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
