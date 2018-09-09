
lineIntersect = require('@turf/line-intersect').default
lineDistance = require('@turf/point-to-line-distance').default
distance = require('@turf/distance').default
centroid = require('@turf/centroid').default
agentmaps = require('agentmaps').default
routing = require('agentmaps/src/routing')

////// dario code ////

map_data.features = map_data.features.filter((item) => item.properties.tags.highway && ['trunk', 'primary', 'secondary', 'tertiary', 'unclassified'].indexOf(item.properties.tags.highway) !== -1)
map = L.map("demo_map").fitBounds(bounding_points);
// map.setZoom(map.getZoom()+1)
L.tileLayer(
    "http://{s}.tiles.wmflabs.org/bw-mapnik/{z}/{x}/{y}.png",
    // "http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
    {
        attribution: "Thanks to <a href=\"http://openstreetmap.org\">OpenStreetMap</a> community",
    }
).addTo(map);

let street_options = {
  "color": "yellow",
  "weight": 4,
  "opacity": .5
};
let unit_options = {
  front_buffer: 6,
  side_buffer: 3,
  length: 10,
  depth: 18,
  color: "blue"
};
// let embassies_coords = embassies.map((item) => {
//   return [item.geometry.y, item.geometry.x]
// });
// embassies_coords = embassies_coords.filter((item) => item[0] > bounding_points[0][0] && item[0] < bounding_points[1][0] && item[1] > bounding_points[1][1] && item[1] < bounding_points[1][1])


agentmap = L.A.agentmap(map);

let options_set = {
  school: {
    "color": "#e34a33",
    "weight": 1,
    "opacity": .87,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 100,
    "depth": 100
  },
  embassy: {
    "color": "green",
    "weight": 1,
    "opacity": .87,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 100,
    "depth": 100
  },
  standard: {
    "color": "grey",
    "weight": 1,
    "opacity": .3,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 25,
    "depth": 25
  },
  person: {
    radius: .5,
    color: "red",
    fillColor: "red",
    opacity: 0.7
  },
  student: {
    radius: .5,
    color: "green",
    fillColor: "green",
    opacity: 0.7
  },
  embassyVip: {
    radius: .5,
    color: "blue",
    fillColor: "blue",
    opacity: 0.7
  }
};

unit_options = options_set.standard
congregate_idxs = [0,1,2]

agentmap.buildingify(bounding_points, map_data, street_options, unit_options, units_data);
congregate_idxs.forEach(function(idx){
    unit = agentmap.units.getLayers()[idx]
    unit.feature.properties.isCongregationPoint = true
    circ = new L.circle(unit.getCenter(), 1, {
      fillOpacity: 0.0,
      opacity: 0.0
    })
    .bindTooltip('EVAC', { permanent: true, direction: "center", offset: [0,10], className: "flagstyle"})
    .addTo(map);
  }
)
agentmap.units.eachLayer(function(unit){
  if(false){//unit.feature.properties.isCongregationPoint){
    //
  } else {
    unittype = unit.feature.properties.poi_type
    unitmarker = false
    if(unittype === "embassy"){
      unitmarker = '🤝🏾'
    } else if (unittype === "school") {
      unitmarker = '🎓'
    }
    if(unitmarker){
      circ = new L.circle(unit.getCenter(), 1, {
        fillOpacity: 0.0,
        opacity: 0.0
      })
      .bindTooltip(unitmarker, { permanent: true, direction: "center", offset: [-10,0], className: "labelstyle"})
      .addTo(map);
    }
  }
})



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

  sel = Math.random()
  if(sel < 0.2){
    poitype = "embassy"
  } else if (sel < 0.4) {
    poitype = "school"
  } else {
    poitype = "standard"
  }

  unitset = this.units.getLayers().filter(unit => unit.feature.properties.poi_type === poitype)

  // let unit = this.units.getLayers()[Math.floor(agentmap.units.count() * Math.random())],
  let unit = unitset[Math.floor(unitset.length * Math.random())],
	unit_id = this.units.getLayerId(unit),
	center_point = centroid(unit.feature);
	center_point.properties.place = {"type": "unit", "id": unit_id},
	center_point.properties.layer_options = {radius: .5, color: "red", fillColor: "red"};

  if (unit.feature.properties.poi_type === 'school') {
    center_point.properties.layer_options = Math.random() > 0.5 ? options_set.person : options_set.student;
  } else if (unit.feature.properties.poi_type === 'embassy') {
    center_point.properties.layer_options = Math.random() > 0.8 ? options_set.person : options_set.embassyVip;
  } else {
    center_point.properties.layer_options = options_set.person
  }

  return center_point;

}


window.numAgents = 500
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
                agent.setTravelToPlace(closest_unit_center, {type: "unit", id: agentmap.units.getLayerId(closest_unit)}, 4, true, true);
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
