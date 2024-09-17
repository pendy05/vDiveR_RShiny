// disable the tab second host
shinyjs.disableTab = function (name) {
  var tab = $(".nav-tabs li a[data-value=" + name + "]");
  tab.bind("click.tab", function (e) {
    e.preventDefault();
    return false;
  });
  tab.addClass("disabled");
};

// enable the tab second host
shinyjs.enableTab = function (name) {
  var tab = $(".nav-tabs li a[data-value=" + name + "]");
  tab.unbind("click.tab");
  tab.removeClass("disabled");
};
