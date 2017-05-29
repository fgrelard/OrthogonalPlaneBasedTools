function changeImage(element) {
    var pathOurs = "Ours/"+element;
    var pathCornea = "Cornea/"+element;
    var pathLee = "Lee/"+element;
    var pathPalagyi = "Palagyi/"+element;
    var pathCouprie = "Couprie/"+element;
    document.getElementById('imgOurMethod').src = "Ours/"+element;
    document.getElementById('largeOurs').style.backgroundImage = "url("+pathOurs+")";
    document.getElementById('largeOurs').style.backgroundRepeat = "no-repeat";
    document.getElementById('cornea').src = "Cornea/"+element;
    document.getElementById('largeCornea').style.backgroundImage = "url("+pathCornea+")";
    document.getElementById('largeCornea').style.backgroundRepeat = "no-repeat";
    document.getElementById('lee').src = "Lee/"+element;
    document.getElementById('largeLee').style.backgroundImage = "url("+pathLee+")";
    document.getElementById('largeLee').style.backgroundRepeat = "no-repeat";
    document.getElementById('palagyi').src = "Palagyi/"+element;
    document.getElementById('largePalagyi').style.backgroundImage = "url("+pathPalagyi+")";
    document.getElementById('largePalagyi').style.backgroundRepeat = "no-repeat";
    document.getElementById('couprie').src = "Couprie/"+element;
    document.getElementById('largeCouprie').style.backgroundImage = "url("+pathCouprie+")";
    document.getElementById('largeCouprie').style.backgroundRepeat = "no-repeat";
}
