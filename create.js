let map = L.map("demo_map").fitBounds(bounding_points);
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
	center_point = turf.centroid(unit.feature);
	center_point.properties.place = {"type": "unit", "id": unit_id},
	center_point.properties.layer_options = {radius: .5, color: "red", fillColor: "red"};
  return center_point;

}

reset_from_saved = false
function restoreFromSaveAgentMaker(id){
	let index = this.agents.count();
  center_point = saved_points[index]
	return center_point;
}

function saveAgents(){
  saved_points = []
  saved_goals = []
  agentmap.agents.eachLayer(function(a){
    m = L.marker(a._latlng)
    center_point = m.toGeoJSON()
    center_point.properties.place = a.place
    center_point.properties.layer_options = a.options
    saved_points.push(center_point)
    if (a.trip.path.length > 0){
      goal = a.trip.path[a.trip.path.length-1].new_place
    } else {
      goal = null
    }
    saved_goals.push(goal)
  })
}

let numAgents = 20
idxarr = [...Array(numAgents).keys()];
shuffleArray(idxarr)
agentmap.agentify(numAgents, randomUnitAgentMaker);
// agentmap.agentify(numAgents, agentmap.seqUnitAgentMaker);

agentmap.controller = function() {
    if (reset_from_saved) {
      // can we figure out how to match the saved goals with the new agents?
      reset_from_saved = false
    }
    if (agentmap.state.ticks === 0) {
        agentmap.agents.eachLayer(function(agent) {
            if (agent.trip.path.length == 0){
              // let random_index = Math.floor(agentmap.units.count() * Math.random()),
              let random_index = 0,
              random_unit = agentmap.units.getLayers()[random_index],
              random_unit_id = agentmap.units.getLayerId(random_unit),
              random_unit_center = random_unit.getBounds().getCenter();

              agent.scheduleTrip(random_unit_center, {type: "unit", id: random_unit_id});
            }
        })
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


recreateAgentmap = function(){
  saveAgents()
  agentmap.clear()
  agentmap = L.A.agentmap(map);
  agentmap.buildingify(bounding_points, map_data, street_options, unit_options, units_data);
  agentmap.agentify(numAgents, restoreFromSaveAgentMaker);
  // agentmap.agentify(numAgents, agentmap.seqUnitAgentMaker);

  agentmap.controller = function() {
      if (reset_from_saved) {
        // can we figure out how to match the saved goals with the new agents?
        reset_from_saved = false
      }
      if (agentmap.state.ticks === 0) {
          agentmap.agents.eachLayer(function(agent) {
              if (agent.trip.path.length == 0){
                // let random_index = Math.floor(agentmap.units.count() * Math.random()),
                let random_index = 0,
                random_unit = agentmap.units.getLayers()[random_index],
                random_unit_id = agentmap.units.getLayerId(random_unit),
                random_unit_center = random_unit.getBounds().getCenter();

                agent.scheduleTrip(random_unit_center, {type: "unit", id: random_unit_id});
              }
          })
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
}
deleteStreet = function(street){
  map_data.features = map_data.features.filter(function(feature){feature.id != street.feature.id})
}
handleClick = function(event){
  deleteStreet(event.target)
  recreateAgentmap()
}
// agentmap.streets.eachLayer(function(unit){unit.once("click", function(street){agentmap.streets.removeLayer(street)})})
agentmap.streets.eachLayer(function(street){street.on("click", handleClick)})
