(defproject com.quantisan/clatrix "0.2.2"
  :description "Because using matrices in Clojure needs to not suck."
  :url "http://quantisan.github.com/clatrix"
  :license {:name "MIT License"
            :url "http://www.opensource.org/licenses/MIT"}
  :resource-paths ["native"]
  :plugins [[lein-expectations "0.0.8"]]
  :dev-dependencies [[lein-expectations "0.0.8"]
                     [expectations "1.4.16"]]
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [slingshot "0.8.0"]
                 [org.jblas/jblas "1.2.2"]]
  :profiles {:dev {:dependencies [[expectations "1.4.16"]]}})